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

var km = require('../js/index');

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

function rRes(str) {
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
function runR(tte, ev) {
	var file = tmp.fileSync();
	fs.writeSync(file.fd, rSurvScript(tte, ev));
	return rRes(cp.execSync(`R --slave < ${file.name}`, {
		encoding: 'utf8',
		stdio: ['pipe', 'pipe', 'ignore']
	}));
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

describe('kaplan-meier', function () {
	describe('#compute', function () {

		property('should match survival', 'array (integer & bool)', function (arr) {
			if (arr.length > 0) {
				var tte = _.pluck(arr, 0);
				var ev = _.pluck(arr, 1);
				// Why is rate in km.compute?
				var res = km.compute(tte, ev).map(r => _.omit(r, ['rate', 'e']));
				var Rres = runR(tte, ev);
			}
			// Using an assert instead of returning a boolean, so we can
			// report the invariant that failed.
			assert(kmEqual(res, Rres, precision),
				`js ${JSON.stringify(res)} deepEqual R ${JSON.stringify(Rres)}`);
			return true;
		});
	});
});
