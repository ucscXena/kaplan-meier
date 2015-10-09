/*global describe: false, require: false, it: false, console: false, child_process: false */
'use strict';
var _ = require('underscore');
var fs = require('fs');
var assert = require('assert');
var tmp = require('tmp');
var jsc = require('jsverify');
var jscOpts = {
	size: 50,
	tests: 100
};
var cp = require('child_process');

var km = require('../js/index').init(_);

function property(name, prop, fn) {
	it(name, function () {
		jsc.assert(jsc.forall(prop, fn), jscOpts);
	});
}

function rSurvScript(tte, ev) {
	var evb = ev.map(b => b.toString().toUpperCase());
	return [
		'library(survival)',
		`time <- c(${tte})`,
		`status <- c(${evb})`,
		'a <- survfit(Surv(time, status)~1)',
		"write.table(data.frame(a$n.risk, a$n.event, a$time, a$surv), file=stdout(), sep='\\t')"
	].join('\n');
}

function chomp(str) {
	return str[str.length - 1] === '\n' ? str.slice(0, str.length - 1) : str;
}

function rSurvResp(str) {
	return _.map(
			chomp(str).split('\n').slice(1)      // split rows & drop the header
				.map(l => l.split('\t').slice(1) // split columns & drop the row index
						.map(s => parseFloat(s, 10))),
			r => _.object(['n', 'd', 't', 's'], r));
}

// The best we can do for now is write R scripts to files &
// execute them. This is very slow. The C API looks like a
// rat hole. Perhaps try this:
//
// http://stackoverflow.com/questions/9370609/piping-stdin-to-r
function runR(script) {
	var file = tmp.fileSync();
	fs.writeSync(file.fd, script);
	return cp.execSync(`R --slave < ${file.name}`, {
		encoding: 'utf8',
		stdio: ['pipe', 'pipe', 'ignore']
	});
}

function cmpPt(a, b, p) {
	return a.n === b.n &&
		a.d === b.d &&
		a.t === b.t &&
		Math.abs(a.s - b.s) < p;
}

function kmEqual(a, b, p) {
	return _.every(a, (ai, i) => cmpPt(ai, b[i], p));
}

var precision = 0.00001;

var subjectsToTimes = subjects => ({
	tte: _.pluck(subjects, 0),
	ev: _.pluck(subjects, 1)
});

var repeat = (i, n) => _.range(n).map(() => i);

// groups: [{tte, ev}, ...]
function rLogRankScript(groups) {
	var tagged = groups.map(({tte, ev}, i) =>
			({tte: tte, ev: ev, group: repeat(i, tte.length)})),
		tte = _.flatten(_.pluck(tagged, 'tte')),
		ev = _.flatten(_.pluck(tagged, 'ev')).map(b => b.toString().toUpperCase()),
		group = _.flatten(_.pluck(tagged, 'group'));
	return [
		'library(survival)',
		`time <- c(${tte})`,
		`status <- c(${ev})`,
		`group <- c(${group})`,
		'x <- survdiff(Surv(time, status)~group)',
		'if (length(x$n)==1)  {',
		'	p <- 1-pchisq(x$chisq, 1)',
		'} else {',
		'	if (is.matrix(x$obs)){',
		'		etmp <- apply(x$exp,1,sum)',
		'	} else {',
		'		etmp <- x$exp',
		'	}',
		'	df <- (sum(1*(etmp>0))) -1',
		'	p <- 1-pchisq(x$chisq, df)',
		'}',
		'write(p, file=stdout())'
	].join('\n');
}

describe('kaplan-meier', function () {
	describe('#compute', function () {
		property('should match R.survival', 'array (integer & bool)', function (arr) {
			var tte, ev, res, Rres;
			if (arr.length > 0) {
				tte = _.pluck(arr, 0);
				ev = _.pluck(arr, 1);
				// Why is rate in km.compute?
				res = km.compute(tte, ev).map(r => _.omit(r, ['rate', 'e']));
				Rres = rSurvResp(runR(rSurvScript(tte, ev)));
			}
			// Using an assert instead of returning a boolean, so we can
			// report the invariant that failed.
			assert(kmEqual(res, Rres, precision),
				`js ${JSON.stringify(res)} deepEqual R ${JSON.stringify(Rres)}`);
			return true;
		});
	});

	describe('#logranktest', function () {
		property('should match R.survival',
			// This construction ensures we have at least 2 groups. nearray
			// enforces that the groups are not empty.
			'nearray (integer & bool) & nearray(nearray (integer & bool))',
			function ([group0, groupsN]) {
				var groups = [group0].concat(groupsN),
					groupTable = groups.map(subjectsToTimes),
					logrank = km.logranktest(groupTable),
					Rlogrank = NaN;

				try {
					Rlogrank = parseFloat(runR(rLogRankScript(groupTable)));
				} catch (e) {
					if (!isNaN(logrank)) { // not expected
						throw e;           // re-throw
					}
				}

				// Using an assert instead of returning a boolean, so we can
				// report the invariant that failed.
				assert(Math.abs(logrank.pValue - Rlogrank) < precision ||
					isNaN(Rlogrank) && isNaN(logrank),
					`| ${logrank.pValue} (js) - ${Rlogrank} (R) | < ${precision}`);
				return true;
		});
	});
});
