/*global require: false, module: false */
'use strict';

var config = require('./webpack.config');

config.output.filename = "testBundle.js";
config.entry = './test/all.js';
module.exports = config;
