import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

// GitHub project pages are served below /<repository-name>/. The deployment
// workflow supplies that path; local development and custom-domain builds use /.
export default defineConfig({
  base: process.env.VITE_BASE_PATH || '/',
  plugins: [react()],
});

