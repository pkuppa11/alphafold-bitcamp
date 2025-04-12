<script>
    export let top_hit; // Passed from your router/navigation
    
    let isLoading = true;
    let embedUrl = `https://alphafold.ebi.ac.uk/entry/${top_hit}`;
    
    // Optional: Fetch additional metadata
    let metadata = null;
    
    onMount(async () => {
      try {
        // You could fetch additional protein info here if needed
        // const response = await fetch(`/api/protein-info?accession=${top_hit}`);
        // metadata = await response.json();
      } finally {
        isLoading = false;
      }
    });
  </script>
  
  <style>
    .model-container {
      width: 100%;
      height: 80vh;
      display: flex;
      flex-direction: column;
      gap: 1rem;
    }
    
    iframe {
      border: none;
      border-radius: 8px;
      box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
      flex-grow: 1;
    }
    
    .header {
      display: flex;
      justify-content: space-between;
      align-items: center;
    }
    
    .loading {
      display: flex;
      justify-content: center;
      align-items: center;
      height: 100%;
    }
    
    .error {
      color: #ef4444;
      text-align: center;
    }
  </style>
  
  <div class="model-container">
    {#if isLoading}
      <div class="loading">
        <p>Loading AlphaFold model for {top_hit}...</p>
      </div>
    {:else}
      <div class="header">
        <h2>AlphaFold Prediction: {top_hit}</h2>
        <a href={embedUrl} target="_blank" rel="noopener">
          Open in Full Screen →
        </a>
      </div>
      
      <iframe
        title={`AlphaFold structure for ${top_hit}`}
        src={embedUrl}
        width="100%"
        height="100%"
        allowfullscreen
      ></iframe>
      
      {#if metadata}
        <div class="metadata">
          <!-- Display any additional protein info here -->
        </div>
      {/if}
    {/if}
  </div>