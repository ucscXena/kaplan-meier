/*global require: false, module: false */
'use strict';

var _ = require('underscore');

var jStat = require('jStat').jStat;

var reduce = _.reduce,
	map = _.map,
	groupBy = _.groupBy,
	sortBy = _.sortBy,
	last = _.last,
	uniq = _.uniq,
	pluck = _.pluck,
	filter = _.filter;

// jStat methods annoyingly demote types if they have length one. This makes
// them fail to compose with other methods.  Here we re-assert the proper types
// for multiply and transpose.
function multiply(a, b) {
	var r = jStat.multiply(a, b);
	return r.length ? r : [[r]];
}

function transpose(a) {
	var r = jStat.transpose(a);
	return r[0].length ? r : [r];
}

// Patch for bug, in case we need it. Affects reuse of matrices that
// are passed to inv, aug, gauss_jordan.
// https://github.com/jstat/jstat/issues/167
//jStat.aug = (a, b) => a.map((row, i) => row.concat(b[i]));

// Compute at-risk, exiting, and deaths for each time t_i, from
// a list of events.
// tte: [number, ...]
// ev:  [boolean, ...]
// returns: [{n, e, d, t}, ...]
function timeTable(tte, ev) {
	var exits = sortBy(map(tte, (x, i) => ({tte: x, ev: ev[i]})), 'tte'), // sort and collate
		uexits = uniq(pluck(exits, 'tte'), true),                // unique tte
		gexits = groupBy(exits, x => x.tte);                     // group by common time of exit
	return reduce(uexits, function (a, tte) {                // compute d_i, n_i for times t_i (including censor times)
		var group = gexits[tte],
		l = last(a) || {n: exits.length, e: 0},
		events = filter(group, x => x.ev);

		a.push({
			n: l.n - l.e,     // at risk
			e: group.length,  // number exiting
			d: events.length, // number events (death)
			t: group[0].tte   // time
		});
		return a;
	}, []);
}

// kaplan-meier
// See http://en.wikipedia.org/wiki/Kaplan%E2%80%93Meier_estimator
//
// tte  time to exit (event or censor)
// ev   is truthy if there is an event.
function compute(tte, ev) {
	var dini = timeTable(tte, ev),

		// s : the survival probability from t=0 to the particular time (i.e. the end of the time interval)
		// rate : the chance of an event happened within the time interval (as in t and the previous t with an event)
		si = reduce(dini, function (a, dn) { // survival at each t_i (including censor times)
			var l = last(a) || { s: 1 };
			if (dn.d) {                      // there were events at this t_i
				a.push({t: dn.t, e: true, s: l.s * (1 - dn.d / dn.n), n: dn.n, d: dn.d, rate: dn.d / dn.n});
			} else {                         // only censors
				a.push({t: dn.t, e: false, s: l.s, n: dn.n, d: dn.d, rate: null});
			}
			return a;
		}, []);

	return si;
}


//log-rank test of the difference between KM plots

// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059453/
// a good article to understand KM and comparing KM plots using log-rank test,
// they used the pearson chisquared test to compute test statistics
// sum of (O-E)^2/E

// http://oto.sagepub.com/content/143/3/331.long
// a good article to understand KM and comparing KM plots using log-rank test and hazardous ratio test
// they also used the pearson chisquared test to compute test statistics

// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC403858/
// introduce pearson chi-square to compute logrank statistics, however mentioned there is another way

// http://ssp.unl.edu/Log%20Rank%20Test%20For%20More%20Than%202%20Groups.pdf
// gives basic idea of the "other" way R seems to use the "other" way
// (O-E)^2/V V is variance for two groups and covariance for multiple groups

// https://cran.r-project.org/web/packages/survival/survival.pdf
// R use (O-E)^2/V V is variance for two groups and covariance for multiple groups

//https://github.com/CamDavidsonPilon/lifelines/blob/master/lifelines/statistics.py
//python implementation, identical results to R

// covariance calculation
// https://books.google.com/books?id=nPkjIEVY-CsC&pg=PA451&lpg=PA451&dq=multivariate+hypergeometric+distribution+covariance&source=bl&ots=yoieGfA4bu&sig=dhRcSYKcYiqLXBPZWOaqzciViMs&hl=en&sa=X&ved=0CEQQ6AEwBmoVChMIkqbU09SuyAIVgimICh0J3w1x#v=onepage&q=multivariate%20hypergeometric%20distribution%20covariance&f=false

//https://plot.ly/ipython-notebooks/survival-analysis-r-vs-python/#Using-R
// R online tutorial

// chisquare distribution at
// https://github.com/jstat/jstat/blob/master/src/distribution.js
// testing jStat accuracy: http://www.socscistatistics.com/pvalues/chidistribution.aspx

// p value = 1- jStat.chisquare.cdf(x, dof );  -- x is chisquare statistics, dof is degree of freedom
// for comparing two plots, the dof is n-1 = 1, comparing three plots dof = n-1 = 2

// given a theoretical survival curve (si), and tte + ev ( tte and ev is the data ),
// compute the expected total number of events
// report observed n events, expected n events. pearson's chi-square component (O-E)^2/E

function expectedObservedEventNumber(si, tte, ev) {
	var data = timeTable(tte, ev),
		expectedNumber,
		observedNumber,
		dataByTimeTable = [];

	si = si.filter(item => item.e);

	expectedNumber = reduce(si, function (memo, item) {
		var pointerInData = _.find(data, x => x.t >= item.t);

		if (pointerInData) {
			var expected = pointerInData.n * item.rate;
			dataByTimeTable.push(pointerInData);
			return memo + expected;
		}
		else {
			return memo;
		}

	}, 0);

	observedNumber = filter(ev, x => x).length;

	return {
		expected: expectedNumber,
		observed: observedNumber,
		dataByTimeTable: dataByTimeTable,
		timeNumber: dataByTimeTable.length
	};
}

function covariance(allGroupsRes, OETable) {
	var vv = jStat.zeros(OETable.length),
		i, j, //groups
		t, //timeIndex
		N, //total number of samples
		Ki, Kj, // at risk number from each group
		n; //total observed

	for (i = 0; i < OETable.length; i++) {
		for (j = i; j < OETable.length; j++) {
			for (t = 0; t < allGroupsRes.length; t++) {
				N = allGroupsRes[t].n;
				n = allGroupsRes[t].d;
				if (t < OETable[i].timeNumber && t < OETable[j].timeNumber) {
					Ki = OETable[i].dataByTimeTable[t].n;
					Kj = OETable[j].dataByTimeTable[t].n;
					// https://books.google.com/books?id=nPkjIEVY-CsC&pg=PA451&lpg=PA451&dq=multivariate+hypergeometric+distribution+covariance&source=bl&ots=yoieGfA4bu&sig=dhRcSYKcYiqLXBPZWOaqzciViMs&hl=en&sa=X&ved=0CEQQ6AEwBmoVChMIkqbU09SuyAIVgimICh0J3w1x#v=onepage&q=multivariate%20hypergeometric%20distribution%20covariance&f=false
					// when N==1: only 1 subject, no variance
					if (i !== j && N !== 1) {
						vv[i][j] -= n * Ki * Kj * (N - n) / (N * N * (N - 1));
						vv[j][i] = vv[i][j];
					}
					else if (N !== 1) {  // i==j
						vv[i][i] += n * Ki * (N - Ki) * (N - n) / (N * N * (N - 1));
					}
				}
			}
		}
	}
	return vv;
}

// allGroupsRes: km of all groups combined?
// groupedDataTable: [{tte, ev}, ...]
function logranktest (allGroupsRes, groupedDataTable) {
	var KMStats,
		pValue,
		dof, // degree of freedom
		OETable,
		OMinusEVector, OMinusEVectorMinus1,// O-E and O-E drop the last element
		vv, vvMinus1; //covariant matrix and covraiance matrix drops the last row and column

	OETable = groupedDataTable
				.map(({tte, ev}) => expectedObservedEventNumber(allGroupsRes, tte, ev))
				.filter(r => r.expected);

	OMinusEVector = map(OETable, r => r.observed - r.expected);
	// logrank stats covariance matrix vv
	vv = covariance(allGroupsRes, OETable);

	OMinusEVectorMinus1 = OMinusEVector.slice(1); // Drop one dimension from O-E and variance
	vvMinus1 = vv.slice(1).map(r => r.slice(1));

	dof = OETable.length - 1;

	var chi = 0;
	if (dof > 0) {
		var m = [OMinusEVectorMinus1],
			mT = transpose(m),
			vvMinus1Inv = jStat.inv(vvMinus1),
			mfinal = multiply(multiply(m, vvMinus1Inv), mT);

		KMStats = mfinal[0][0];

		chi = jStat.chisquare.cdf(KMStats, dof);
	}

	pValue = 1 - chi;
	return {
		dof: dof,
		KMStats: KMStats,
		pValue: pValue
	};
}

module.exports = {
	compute: compute,
	expectedObservedEventNumber: expectedObservedEventNumber,
	logranktest: logranktest
};
