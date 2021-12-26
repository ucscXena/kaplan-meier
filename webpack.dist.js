var config = require('./webpack.config');

config.output = {
    path: 'dist',
    filename: 'kaplan-meier.js',
    library: 'kaplan-meier',
    libraryTarget: 'umd'
};

config.externals = [
	'jstat'
];

module.exports = config;
