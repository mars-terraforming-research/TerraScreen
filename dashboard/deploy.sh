rm -rf dist docs
npm ci    # or: npm install
npm run build                                             
git add -A
git commit -m "Build to docs with relative asset URLs"
git push

