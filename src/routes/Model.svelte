<script>
  import { page } from '$app/stores';
  
  let isLoading = true;
  let error = null;
  let modelUrl = '';

  $: {
    const accession = $page.url.searchParams.get('accession');
    if (!accession) {
      error = "No protein accession provided";
    } else {
      modelUrl = `https://alphafold.ebi.ac.uk/entry/${accession}`;
    }
    isLoading = false;
  }
</script>

<div class="model-container">
  {#if isLoading}
    <div class="loading">Loading structure...</div>
  
  {:else if error}
    <div class="error">
      {error}
      <a href="/" class="back-link">← Back to Predictor</a>
    </div>
  
  {:else}
    <div class="header">
      <h1>AlphaFold Structure: {$page.url.searchParams.get('accession')}</h1>
      <a href="/" class="back-button">← Back</a>
    </div>
    
    <div class="viewer-container">
      <iframe
        title="AlphaFold Structure Viewer"
        src={modelUrl}
        allowfullscreen
      ></iframe>
    </div>
  {/if}
</div>

<style>
  .model-container {
    width: 100%;
    height: 100vh;
    display: flex;
    flex-direction: column;
    background: linear-gradient(135deg, #0f172a, #1e3a8a);
    color: white;
    font-family: 'Roboto Mono', monospace;
  }

  .header {
    padding: 1.5rem;
    background: rgba(255, 255, 255, 0.1);
    display: flex;
    justify-content: space-between;
    align-items: center;
  }

  h1 {
    font-family: 'Orbitron', sans-serif;
    margin: 0;
    background: linear-gradient(to right, #38bdf8, #4ade80);
    -webkit-background-clip: text;
    background-clip: text;
    color: transparent;
    font-size: 1.5rem;
  }

  .viewer-container {
    flex-grow: 1;
    position: relative;
  }

  iframe {
    width: 100%;
    height: 100%;
    border: none;
  }

  .loading, .error {
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
    height: 100%;
    gap: 1rem;
  }

  .error {
    color: #ff6b6b;
  }

  .back-button, .back-link {
    background: rgba(255, 255, 255, 0.1);
    border: 1px solid rgba(255, 255, 255, 0.3);
    color: white;
    padding: 0.5rem 1rem;
    border-radius: 50px;
    text-decoration: none;
    transition: all 0.3s ease;
  }

  .back-button:hover, .back-link:hover {
    background: rgba(255, 255, 255, 0.2);
  }
</style>