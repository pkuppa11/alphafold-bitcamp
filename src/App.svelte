<script>
  import { Router, Route } from 'svelte-routing';
  import Home from './routes/Home.svelte';
  import Model from './routes/Model.svelte';
  import { writable } from 'svelte/store';

  const currentUrl = writable(window.location.hash.slice(1) || '/');

  window.addEventListener('hashchange', () => {
    currentUrl.set(window.location.hash.slice(1) || '/');
  });
</script>

<Router url = {$currentUrl} >
  <Route path="/" component={Home} />
  <Route path="/model" component={Model} />
</Router>