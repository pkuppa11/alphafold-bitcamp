<script>
  import { onMount } from "svelte";

  let accession = "";
  let isLoading = true;
  let error = "";

  onMount(() => {
    try {
      const hash = window.location.hash;
      const queryString = hash.split("?")[1];
      const params = new URLSearchParams(queryString);
      accession = params.get("accession");

      if (!accession) {
        error = "No accession code found in URL.";
      }
    } catch (e) {
      error = "Failed to parse URL.";
    } finally {
      isLoading = false;
    }
  });
</script>

<main>
  <div class="content">
    {#if isLoading}
      <p class="status">🔄 Loading protein structure...</p>
    {:else if error}
      <p class="status error">⚠️ {error}</p>
    {:else}
      <h2>Protein Structure: <span class="accession">{accession}</span></h2>
      <div class="viewer-wrapper">
        <iframe
          title="AlphaFold Structure Viewer"
          src={`https://alphafold.ebi.ac.uk/entry/${accession}`}
          width="100%"
          height="600px"
          allowfullscreen
        />
      </div>
    {/if}
  </div>
</main>

<style>
  :global(body) {
    margin: 0;
    font-family: "Roboto Mono", monospace;
    background: linear-gradient(135deg, #0f172a, #1e3a8a);
    color: white;
  }

  main {
    min-height: 100vh;
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 3rem 2rem;
    position: relative;
    overflow: hidden;
  }

  .content {
    max-width: 60rem;
    width: 100%;
    text-align: center;
    z-index: 2;
    background: rgba(255, 255, 255, 0.05);
    padding: 2rem;
    border-radius: 1rem;
    box-shadow: 0 10px 40px rgba(0, 0, 0, 0.3);
    backdrop-filter: blur(12px);
    border: 1px solid rgba(255, 255, 255, 0.1);
  }

  h2 {
    font-family: "Orbitron", sans-serif;
    font-size: 2rem;
    margin-bottom: 1.5rem;
    letter-spacing: 1px;
    background: linear-gradient(to right, #38bdf8, #4ade80);
    -webkit-background-clip: text;
    background-clip: text;
    color: transparent;
  }

  .accession {
    font-weight: 600;
    color: #a5f3fc;
    text-shadow: 0 0 10px rgba(56, 189, 248, 0.6);
  }

  .status {
    font-size: 1.2rem;
    padding: 1rem;
    color: #cbd5e1;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 12px;
    display: inline-block;
    animation: fade-in-up 0.6s ease-out forwards;
  }

  .status.error {
    color: #fca5a5;
    background: rgba(239, 68, 68, 0.1);
  }

  .viewer-wrapper {
    border-radius: 12px;
    overflow: hidden;
    box-shadow: 0 0 25px rgba(56, 189, 248, 0.3);
    animation: fade-in-up 0.8s 0.2s ease-out forwards;
  }

  iframe {
    border: none;
    width: 100%;
    height: 600px;
  }

  @keyframes fade-in-up {
    from {
      opacity: 0;
      transform: translateY(30px);
    }
    to {
      opacity: 1;
      transform: translateY(0);
    }
  }
</style>
