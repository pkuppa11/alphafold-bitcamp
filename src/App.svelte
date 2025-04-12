<script>
  let sequence = "";
</script>

<style>
  main {
    min-height: 100vh;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    padding: 2rem;
    background: linear-gradient(to bottom right, #3b82f6, #10b981);
    color: white;
    font-family: sans-serif;
  }

  h1 {
    font-size: 3rem;
    font-weight: 800;
    margin-bottom: 2rem;
    animation: fade-in-up 1s ease-out forwards;
    text-shadow: 2px 2px 8px rgba(0, 0, 0, 0.3);
  }

  input {
    width: 100%;
    max-width: 40rem;
    padding: 1rem 2rem;
    border-radius: 9999px;
    font-size: 1.1rem;
    color: black;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
    border: none;
    outline: none;
    transition: box-shadow 0.3s ease, transform 0.3s ease;
    animation: slide-in 1s ease-out forwards;
  }

  input:focus {
    box-shadow: 0 0 0 4px rgba(16, 185, 129, 0.4);
    transform: scale(1.02);
  }

  @keyframes fade-in-up {
    0% {
      opacity: 0;
      transform: translateY(20px);
    }
    100% {
      opacity: 1;
      transform: translateY(0);
    }
  }

  @keyframes slide-in {
    0% {
      opacity: 0;
      transform: translateX(-30px);
    }
    100% {
      opacity: 1;
      transform: translateX(0);
    }
  }
</style>

<main>
  <h1>Alphafold App</h1>
  <input bind:value={sequence} placeholder="Input protein sequence" />
</main>
