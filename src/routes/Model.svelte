<script>
  import { onMount } from "svelte";

  let accession = "";
  let isLoading = true;
  let error = "";

  onMount(() => {
    try {
      // Example URL: http://localhost:8080/#/model?accession=P24941
      const hash = window.location.hash; // → "#/model?accession=P24941"
      const queryString = hash.split('?')[1]; // → "accession=P24941"

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
  {#if isLoading}
    <p>Loading...</p>
  {:else if error}
    <p style="color: red;">{error}</p>
  {:else}
    <h2>Structure for Accession: {accession}</h2>
    <iframe
      title="AlphaFold Structure Viewer"
      src={`https://alphafold.ebi.ac.uk/entry/${accession}`}
      width="100%"
      height="700px"
      style="border: none; margin-top: 1rem;"
    />
  {/if}
</main>

<style>
  main {
    padding: 2rem;
    font-family: "Roboto Mono", monospace;
    background-color: #f0f4f8;
    min-height: 100vh;
  }

  h2 {
    font-size: 1.5rem;
    color: #1e293b;
  }

  p {
    font-size: 1.2rem;
    color: #334155;
  }
</style>
