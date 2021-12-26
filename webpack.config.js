/*global require: false, module: false, __dirname: false */
'use strict';

const path = require('path');

module.exports = {
	//historyApiFallback: true,
	mode: 'production',
	entry: "./js/index",
	output: {
		path: path.join(__dirname, "build"),
		publicPath: "/",
		filename: "[name].js"
	},
	module: {
		rules: [
			{
				test: /\.js$/,
				exclude: /node_modules/,
				use: {
					loader: 'babel-loader',
					options: {
						presets: ['@babel/preset-env'],
						cacheDirectory: true
					}
				}
			}
		]
	},
	plugins: [],
	resolve: {
		extensions: ['', '.js']
	}
};
