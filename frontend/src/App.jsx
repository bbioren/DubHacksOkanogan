import React, { useState } from 'react'

function App() {
  const [query, setQuery] = useState('')
  const [results, setResults] = useState(null)

  const runGeneBot = async () => {
    try {
      const res = await fetch('http://localhost:8000/run_genebot', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query }),
      })
      if (!res.ok) throw new Error(await res.text())
      const data = await res.json()
      setResults(data)
    } catch (err) {
      setResults({ error: String(err) })
    }
  }

  return (
    <div className="container">
      <h1 className="title">🍏 Apple GeneBot</h1>
      <textarea
        className="input"
        value={query}
        onChange={(e) => setQuery(e.target.value)}
        placeholder="Enter apple tree traits..."
      />
      <button className="btn" onClick={runGeneBot}>
        Run GeneBot
      </button>

      {results && (
        <div className="results">
          {results.error ? (
            <pre className="error">{results.error}</pre>
          ) : (
            <>
              <h2>Traits:</h2>
              <p>{Array.isArray(results.traits) ? results.traits.join(', ') : String(results.traits)}</p>
              <h2>Abstracts:</h2>
              <pre>{results.abstracts}</pre>
              <h2>Gene Suggestions:</h2>
              <pre>{results.summary}</pre>
            </>
          )}
        </div>
      )}
    </div>
  )
}

export default App
