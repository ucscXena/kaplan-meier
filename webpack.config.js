/*global require: false, module: false, __dirname: false */
'use strict';
//var HtmlWebpackPlugin = require('html-webpack-plugin');
//var webpack = require('webpack');

module.exports = {
	historyApiFallback: true,
	entry: "./js/index",
	output: {
		path: "build",
		publicPath: "/",
		filename: "[name].js"
	},
	module: {
		loaders: [
			// es7.objectRestSpread failing???
			{ test: /\.js$/, exclude: /node_modules/, loader: 'babel-loader?optional=runtime,cacheDirectory=true'}
		]
	},
	resolve: {
		extensions: ['', '.js'],
		root: __dirname + "/js"
	}
};
