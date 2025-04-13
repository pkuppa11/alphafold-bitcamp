<script>
  import { Router, Route } from 'svelte-routing';
  import Home from './routes/Home.svelte';
  import Model from './routes/Model.svelte';
</script>

<Router>
  <Route path="/" component={Home} />
  <Route path="/model" component={Model} />
</Router>