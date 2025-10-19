// frontend/src/App.jsx
import { useState } from "react";

function App() {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState(null);

  const runGeneBot = async () => {
    const res = await fetch("http://localhost:8000/run_genebot", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ query }),
    });
    const data = await res.json();
    setResults(data);
  };

  return (
    <div className="container mx-auto p-4">
      <h1 className="text-3xl font-bold">🍏 Apple GeneBot</h1>
      <textarea
        className="border p-2 w-full"
        value={query}
        onChange={(e) => setQuery(e.target.value)}
        placeholder="Enter apple tree traits..."
      />
      <button
        className="bg-green-500 text-white px-4 py-2 mt-2 rounded"
        onClick={runGeneBot}
      >
        Run GeneBot
      </button>

      {results && (
        <div className="mt-4">
          <h2 className="font-bold">Traits:</h2>
          <p>{results.traits.join(", ")}</p>
          <h2 className="font-bold">Abstracts:</h2>
          <pre>{results.abstracts}</pre>
          <h2 className="font-bold">Gene Suggestions:</h2>
          <pre>{results.summary}</pre>
        </div>
      )}
    </div>
  );
}

export default App;
