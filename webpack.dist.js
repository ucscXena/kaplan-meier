const path = require('path');

var config = require('./webpack.config');

config.output = {
    path: path.join(__dirname, 'dist'),
    filename: 'kaplan-meier.js',
    library: 'kaplan-meier',
    libraryTarget: 'umd'
};

config.externals = [
	'jStat'
];

module.exports = config;
