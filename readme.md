# kaplan-meier

Basic Kaplan-Meier methods, ported to javascript.

## Dependencies

We assume you're using an es5 shim for array methods.

## Build
The build is based on npm and webpack.
 * Ensure that git and node are installed
   * On OSX, install brew http://brew.sh/
   * `brew install git`
   * `brew install node`
 * `git clone https://github.com/acthp/kaplan-meier.git`
 * `cd kaplan-meier`
 * `npm install`
 * `npm run build`

### Lint

Use `npm run lint` to run the lint rules. We lint with eslint and babel-eslint.

### Test

We test against R. R must be installed.

* `npm test`

### References
 * http://blog.keithcirkel.co.uk/how-to-use-npm-as-a-build-tool/
 * http://webpack.github.io/
 * http://www.youtube.com/watch?v=VkTCL6Nqm6Y
