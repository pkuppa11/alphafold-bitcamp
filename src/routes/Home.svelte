<script>
  import { onMount } from "svelte";
  import { get_root_for_style, get_store_value } from "svelte/internal";
  let sequence = "";
  let particles = [];

  onMount(() => {
    const aminoAcids = [
      "A",
      "R",
      "N",
      "D",
      "C",
      "E",
      "Q",
      "G",
      "H",
      "I",
      "L",
      "K",
      "M",
      "F",
      "P",
      "S",
      "T",
      "W",
      "Y",
      "V",
    ];
    particles = Array.from({ length: 30 }, () => ({
      char: aminoAcids[Math.floor(Math.random() * aminoAcids.length)],
      x: Math.random() * 100,
      y: Math.random() * 100,
      size: Math.random() * 0.8 + 0.5,
      speed: Math.random() * 0.3 + 0.1,
      opacity: Math.random() * 0.4 + 0.2,
    }));
    const animate = () => {
      particles = particles.map((p) => ({
        ...p,
        y: (p.y + p.speed) % 100,
      }));
      requestAnimationFrame(animate);
    };
    animate();
  });
  async function handleFileUpload(e) {
    const file = e.target.files[0];
    if (!file) return;

    try {
      const text = await file.text();
      let content = text
        .split("\n")
        .slice(1)
        .join("")
        .replace(/\s/g, "")
        .toUpperCase();
      content = content.replace(/[^ARNDCEQGHILKMFPSTWYV]/g, "");

      if (content.length === 0) {
        alert("Error: No valid amino acids found in file!");
        return;
      }

      sequence = content;
    } catch (err) {
      alert("Error reading file!");
    }
  }
  function handleInput(e) {
      let value = e.target.value.toUpperCase();
      value = value.replace(/[^ARNDCEQGHILKMFPSTWYV]/g, "");
      sequence = value;
    }

  async function predictStructure() {
    if (!sequence) {
      alert("Please enter a valid sequence.");
      return;
    }

    try {
      const response = await fetch("http://localhost:8000/predict", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ sequence }),
      });

      const data = await response.json();
      if (!response.ok) throw new Error(data.detail || "Unknown error");

      window.location.hash = `model?accession=${data.accession}`;
      alert(`✅ Structure predicted!\nAccession: ${data.accession}`);
      location.reload()
    } catch (err) {
      alert(`❌ Prediction failed: ${err.message}`);
    }
  }

</script>

<main>
  {#each particles as particle, i (i)}
    <div
      class="particle"
      style="
          left: {particle.x}vw;
          top: {particle.y}vh;
          font-size: {particle.size}rem;
          opacity: {particle.opacity};
        "
    >
      {particle.char}
    </div>
  {/each}

  <div class="content">
    <h1>AlphaMutate</h1>
    <p class="subtitle">
      Enter a protein sequence to generate a 3D structure prediction:
    </p>
    <div class="input-container">
      <input
        bind:value={sequence}
        on:input={handleInput}
        placeholder="Example: ACDEFGHIKLMNPQRSTVWYV"
      />
      <label class="upload-button">
        📁 Upload FASTA
        <input
          type="file"
          accept=".fasta,.txt"
          on:change={(e) => handleFileUpload(e)}
          hidden
        />
      </label>
    </div>
    <button class="button" on:click={predictStructure}>Predict Structure →</button>
  </div>
</main>

<style>
  :global(body) {
    margin: 0;
    font-family: "Arial", sans-serif;
  }

  main {
    min-height: 100vh;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    padding: 2rem;
    background: linear-gradient(135deg, #0f172a, #1e3a8a);
    color: white;
    position: relative;
    overflow: hidden;
  }
  .particle {
    position: absolute;
    color: rgba(74, 222, 128, 0.6);
    font-weight: bold;
    user-select: none;
    pointer-events: none;
    z-index: 1;
  }

  .content {
    position: relative;
    z-index: 2;
    width: 100%;
    max-width: 48rem;
    text-align: center;
  }

  h1 {
    font-family: "Orbitron", sans-serif;
    letter-spacing: 1px;
    font-size: 3rem;
    font-weight: 800;
    margin-bottom: 1rem;
    background: linear-gradient(to right, #38bdf8, #4ade80);
    -webkit-background-clip: text;
    background-clip: text;
    color: transparent;
    animation: fade-in-up 0.8s ease-out forwards;
  }

  .subtitle {
    font-size: 1.1rem;
    opacity: 0.8;
    margin-bottom: 2rem;
    animation: fade-in-up 0.8s 0.2s ease-out forwards;
  }
  .input-container {
    display: flex;
    gap: 0.5rem;
    align-items: center;
    margin-bottom: 2rem;
  }

  input {
    width: 100%;
    padding: 1rem 1.5rem;
    border-radius: 50px;
    font-size: 1rem;
    background: rgba(255, 255, 255, 0.1);
    border: 1px solid rgba(255, 255, 255, 0.3);
    color: white;
    box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
    transition: all 0.3s;
  }

  input:focus {
    outline: none;
    border-color: #4ade80;
    box-shadow: 0 0 0 3px rgba(74, 222, 128, 0.3);
  }

  .button {
    background: linear-gradient(to right, #326ecd, #2e814c);
    padding: 1rem 2.5rem;
    border-radius: 50px;
    font-weight: 600;
    border: none;
    color: white;
    cursor: pointer;
    box-shadow: 0 4px 15px rgba(59, 130, 246, 0.3);
    transform: scale(1);
    transition:
      transform 0.2s ease,
      box-shadow 0.2s ease;
  }

  .button:hover {
    transform: scale(1.05);
    box-shadow: 0 6px 20px rgba(59, 130, 246, 0.6);
  }

  .button:active {
    transform: scale(0.98);
    box-shadow: 0 2px 10px rgba(59, 130, 246, 0.4);
  }
  .upload-button {
    background: rgba(255, 255, 255, 0.1);
    border: 1px solid rgba(255, 255, 255, 0.3);
    color: white;
    padding: 0.8rem 1.5rem;
    border-radius: 50px;
    cursor: pointer;
    font-family: "Roboto Mono", monospace;
    font-size: 0.9rem;
    transition: all 0.3s ease;
    white-space: nowrap;
  }

  .upload-button:hover {
    background: rgba(255, 255, 255, 0.2);
  }
  h1 {
    font-family: "Orbitron", sans-serif;
    letter-spacing: 1px;
  }

  .subtitle,
  input,
  .button {
    font-family: "Roboto Mono", monospace;
  }

  @keyframes fade-in-up {
    from {
      opacity: 0;
      transform: translateY(20px);
    }
    to {
      opacity: 1;
      transform: translateY(0);
    }
  }
</style>