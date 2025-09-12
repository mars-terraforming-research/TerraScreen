import { defineConfig } from 'vite';

export default defineConfig({
  base: '/webdash-demo/',   // <-- the repo name, with leading & trailing slashes
  build: { outDir: 'docs' } // optional: write directly to /docs
});
