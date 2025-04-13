<script>
  import { onMount } from "svelte";

  let accessions = [];
  let isLoading = true;
  let error = "";
  let currentMutation = 0;

  onMount(() => {
    try {
      const hash = window.location.hash;
      const queryString = hash.split("?")[1];
      const params = new URLSearchParams(queryString);
      const raw = params.get("accession");

      accessions = raw ? raw.split(",") : [];

      if (accessions.length !== 4) {
        error = "Please provide exactly 4 accession codes separated by commas.";
      }
    } catch (e) {
      error = "Failed to parse accession codes from URL.";
    } finally {
      isLoading = false;
    }
  });

  const prev = () => {
    if (currentMutation > 0) currentMutation -= 1;
  };

  const next = () => {
    if (currentMutation < 2) currentMutation += 1;
  };
</script>

<main>
  <div class="content">
    {#if isLoading}
      <p class="status">🔄 Loading protein structures...</p>
    {:else if error}
      <p class="status error">⚠️ {error}</p>
    {:else}
      <h2>
        🔬 Original Protein: <span class="accession">{accessions[0]}</span>
      </h2>
      <div class="viewer-wrapper main">
        <iframe
          title={`Original Protein - ${accessions[0]}`}
          src={`https://www.ncbi.nlm.nih.gov/Structure/icn3d/?mmdbafid=${accessions[0]}`}
          allowfullscreen
        />
      </div>

      <h3 class="mutations-title">🧪 Mutated Proteins</h3>
      <div class="carousel-with-assistant">
        <div class="carousel">
          <div class="carousel-controls">
            <button
              class="carousel-button"
              on:click={prev}
              disabled={currentMutation === 0}>←</button
            >
            <span
              >Mutation {currentMutation + 1}:
              <span class="accession">{accessions[currentMutation + 1]}</span
              ></span
            >
            <button
              class="carousel-button"
              on:click={next}
              disabled={currentMutation === 2}>→</button
            >
          </div>

          <div class="viewer-wrapper mutation">
            <iframe
              title={`Mutation ${currentMutation + 1} - ${accessions[currentMutation + 1]}`}
              src={`https://www.ncbi.nlm.nih.gov/Structure/icn3d/?mmdbafid=${accessions[currentMutation + 1]}`}
              allowfullscreen
            />
          </div>
          <div class="assistant-box">
            <h4>💬 Protein Assistant</h4>
            <div class="chat-body">
              <p>
                <strong>Assistant:</strong> I'm here to help you analyze mutations!
              </p>
            </div>
          </div>
        </div>
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
}

.content {
  max-width: 90rem;
  width: 100%;
  text-align: center;
  z-index: 2;
  background: rgba(255, 255, 255, 0.05);
  padding: 2.5rem;
  border-radius: 1rem;
  box-shadow: 0 10px 40px rgba(0, 0, 0, 0.3);
  backdrop-filter: blur(12px);
  border: 1px solid rgba(255, 255, 255, 0.1);
}

h2,
h3 {
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
  animation: fade-in-up 0.8s ease-out forwards;
  margin-bottom: 2rem;
}

.viewer-wrapper.main iframe {
  width: 100%;
  height: 600px;
}

.viewer-wrapper.mutation iframe {
  width: 100%;
  height: 500px;
}

.mutations-title {
  font-size: 1.5rem;
  margin: 2rem 0 1rem;
  color: #93c5fd;
}

.carousel-controls {
  display: flex;
  justify-content: center;
  align-items: center;
  gap: 1.5rem;
  margin-bottom: 1rem;
}

.carousel-button {
  padding: 0.6rem 1.2rem;
  font-size: 1.1rem;
  border: none;
  border-radius: 50px;
  background: rgba(255, 255, 255, 0.1);
  color: white;
  cursor: pointer;
  transition: all 0.3s ease;
  font-family: "Roboto Mono", monospace;
}

.carousel-button:hover:not(:disabled) {
  background: rgba(255, 255, 255, 0.25);
}

.carousel-button:disabled {
  opacity: 0.3;
  cursor: not-allowed;
}

.carousel-with-assistant {
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
  gap: 2rem;
  margin-top: 2rem;
  flex-wrap: wrap;
}

.carousel {
  flex: 2 1 60%;
  min-width: 300px;
}

.assistant-box {
  flex: 1 1 35%;
  min-width: 250px;
  background: rgba(255, 255, 255, 0.05);
  border-radius: 12px;
  padding: 1rem;
  border: 1px solid rgba(255, 255, 255, 0.1);
  box-shadow: 0 0 15px rgba(56, 189, 248, 0.2);
  font-family: "Roboto Mono", monospace;
  color: #e0f2fe;
  animation: fade-in-up 0.8s ease-out forwards;
  max-height: 500px;
  overflow-y: auto;
}

.assistant-box h4 {
  font-size: 1.2rem;
  margin-bottom: 1rem;
  background: linear-gradient(to right, #38bdf8, #4ade80);
  -webkit-background-clip: text;
  background-clip: text;
  color: transparent;
}

.chat-body {
  font-size: 0.9rem;
  line-height: 1.4;
  background: rgba(255, 255, 255, 0.03);
  padding: 1rem;
  border-radius: 8px;
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
