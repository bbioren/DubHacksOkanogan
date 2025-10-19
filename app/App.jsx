import React, { useState } from 'react';
import { Search, Loader2, Dna, Microscope, BookOpen } from 'lucide-react';

const AppleGeneBot = () => {
  const [input, setInput] = useState('');
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState(null);
  const [error, setError] = useState('');

  // Preloaded data (same as your Python version)
  const PRELOADED_ABSTRACTS = {
    "wood flexibility": "Genes involved in secondary cell wall biosynthesis: Lignin biosynthesis genes (CCR, CAD), Cellulose Synthase genes (CESA), transcription factors MYB46, NAC. Alterations in these genes modulate wood rigidity and flexibility. CRISPR upregulation increases rigidity, knockout/downregulation increases flexibility.",
    "drooping branches": "Genes affecting branch angle and weeping habit: PIN-FORMED (PIN) auxin transporters, EXPANSINs, TIR1 auxin receptor, LAZY1. Modifying these genes changes branch orientation and drooping.",
    "weeping habit": "Genes regulating branch gravitropism and orientation: LAZY1, TAC1, cytoskeletal regulators. CRISPR edits can induce or reduce weeping forms.",
    "annual fruiting": "Genes controlling flowering and fruiting cycles: FT (Flowering Locus T), TFL1, SOC1, AP1, LFY. CRISPR upregulation promotes yearly fruiting; knockouts can delay flowering.",
    "fruit size": "Genes controlling fruit growth: Cell cycle regulators (CYCD3, CYCB1), Auxin response factors (ARF), Gibberellin pathway genes (GA20ox, GA3ox). Modulating these can increase or reduce fruit size.",
    "sweetness": "Sugar metabolism genes: Sucrose Synthase (SUS), Invertases (INV), Hexose Transporters (HT). Upregulation increases sugar accumulation; knockouts decrease sweetness.",
    "disease resistance": "Resistance genes: R-genes (NBS-LRR family), Pathogenesis-Related proteins (PR1, PR5), WRKY transcription factors. CRISPR can enhance or modify resistance."
  };

  const extractTraits = (text) => {
    // Simple trait extraction (in production, this would call your API)
    const keywords = text.toLowerCase().match(/\b(short|tall|dwarf|sweet|sour|large|small|disease resistant|flexible|drooping|weeping|fruiting|annual|size|sweetness|resistance|wood|branches)\b/g);
    return keywords ? [...new Set(keywords)] : [];
  };

  const getPreloadedAbstracts = (traits) => {
    let abstracts = [];
    traits.forEach(trait => {
      const key = Object.keys(PRELOADED_ABSTRACTS).find(k => 
        k.includes(trait) || trait.includes(k.split(' ')[0])
      );
      if (key) {
        abstracts.push(`**${key.toUpperCase()}**\n${PRELOADED_ABSTRACTS[key]}`);
      }
    });
    return abstracts.length > 0 ? abstracts.join('\n\n---\n\n') : 'No matching abstracts found in preloaded data.';
  };

  const generateGeneSuggestions = (traits, abstracts) => {
    // Simulated gene suggestions (in production, this would call Gemini API)
    const suggestions = [];
    
    if (traits.includes('short') || traits.includes('dwarf')) {
      suggestions.push({
        trait: 'Dwarfing',
        genes: ['GA20ox', 'GA3ox', 'DELLA'],
        function: 'Gibberellin biosynthesis and signaling genes control plant height',
        crispr: 'Knockout GA20ox to reduce gibberellin production → shorter trees'
      });
    }
    
    if (traits.includes('sweet') || traits.includes('sweetness')) {
      suggestions.push({
        trait: 'Sweetness',
        genes: ['SUS1', 'SUS2', 'INV', 'HT1'],
        function: 'Sugar synthesis and transport genes regulate fruit sweetness',
        crispr: 'Upregulate Sucrose Synthase (SUS) genes → increased sugar accumulation'
      });
    }
    
    if (traits.includes('fruiting') || traits.includes('annual')) {
      suggestions.push({
        trait: 'Annual Fruiting',
        genes: ['FT', 'SOC1', 'AP1', 'LFY'],
        function: 'Flowering time genes control reproductive cycles',
        crispr: 'Upregulate FT (Flowering Locus T) → promote consistent yearly fruiting'
      });
    }
    
    if (traits.includes('size') || traits.includes('large')) {
      suggestions.push({
        trait: 'Fruit Size',
        genes: ['CYCD3', 'ARF', 'GA20ox'],
        function: 'Cell cycle and hormone response genes affect fruit growth',
        crispr: 'Upregulate CYCD3 (cell cycle gene) → larger fruit through increased cell division'
      });
    }

    if (suggestions.length === 0) {
      suggestions.push({
        trait: 'General',
        genes: ['Multiple candidates'],
        function: 'Based on your traits, consider consulting specialized literature',
        crispr: 'Multiple CRISPR strategies possible depending on specific goals'
      });
    }

    return suggestions;
  };

  const handleSubmit = async () => {
    if (!input.trim()) {
      setError('Please enter some traits!');
      return;
    }

    setError('');
    setLoading(true);

    // Simulate API delay
    await new Promise(resolve => setTimeout(resolve, 1500));

    try {
      const traits = extractTraits(input);
      const abstracts = getPreloadedAbstracts(traits);
      const geneSuggestions = generateGeneSuggestions(traits, abstracts);

      setResults({
        traits,
        abstracts,
        geneSuggestions
      });
    } catch (err) {
      setError('An error occurred while processing your request.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-green-50 via-white to-green-50">
      <div className="max-w-6xl mx-auto p-6">
        {/* Header */}
        <div className="text-center mb-8">
          <div className="flex items-center justify-center gap-3 mb-4">
            <div className="text-5xl">🍏</div>
            <h1 className="text-4xl font-bold text-green-800">Apple GeneBot</h1>
          </div>
          <p className="text-gray-600 max-w-2xl mx-auto">
            Enter desired apple tree traits in natural language. GeneBot will extract traits, 
            search literature, and suggest candidate genes with potential CRISPR edits.
          </p>
        </div>

        {/* Input Form */}
        <div className="mb-8">
          <div className="bg-white rounded-lg shadow-lg p-6">
            <label className="block text-sm font-medium text-gray-700 mb-2">
              Enter Apple Tree Traits
            </label>
            <textarea
              value={input}
              onChange={(e) => setInput(e.target.value)}
              placeholder="e.g., short trees that fruit annually with sweet, large fruit"
              className="w-full px-4 py-3 border border-gray-300 rounded-lg focus:ring-2 focus:ring-green-500 focus:border-transparent resize-none"
              rows="4"
            />
            <button
              onClick={handleSubmit}
              disabled={loading}
              className="mt-4 w-full bg-green-600 hover:bg-green-700 text-white font-semibold py-3 px-6 rounded-lg transition-colors flex items-center justify-center gap-2 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              {loading ? (
                <>
                  <Loader2 className="animate-spin" size={20} />
                  Processing...
                </>
              ) : (
                <>
                  <Search size={20} />
                  Run GeneBot
                </>
              )}
            </button>
            {error && (
              <p className="mt-2 text-red-600 text-sm">{error}</p>
            )}
          </div>
        </div>

        {/* Results */}
        {results && (
          <div className="space-y-6">
            {/* Extracted Traits */}
            <div className="bg-white rounded-lg shadow-lg p-6">
              <div className="flex items-center gap-2 mb-4">
                <Dna className="text-green-600" size={24} />
                <h2 className="text-xl font-bold text-gray-800">Extracted Traits</h2>
              </div>
              <div className="flex flex-wrap gap-2">
                {results.traits.map((trait, idx) => (
                  <span
                    key={idx}
                    className="px-3 py-1 bg-green-100 text-green-800 rounded-full text-sm font-medium"
                  >
                    {trait}
                  </span>
                ))}
              </div>
            </div>

            {/* Gene Suggestions */}
            <div className="bg-white rounded-lg shadow-lg p-6">
              <div className="flex items-center gap-2 mb-4">
                <Microscope className="text-green-600" size={24} />
                <h2 className="text-xl font-bold text-gray-800">Suggested Genes & CRISPR Edits</h2>
              </div>
              <div className="space-y-4">
                {results.geneSuggestions.map((suggestion, idx) => (
                  <div key={idx} className="border-l-4 border-green-500 pl-4 py-2">
                    <h3 className="font-bold text-gray-800 mb-1">{suggestion.trait}</h3>
                    <p className="text-sm text-gray-600 mb-1">
                      <span className="font-semibold">Genes:</span> {suggestion.genes.join(', ')}
                    </p>
                    <p className="text-sm text-gray-600 mb-1">
                      <span className="font-semibold">Function:</span> {suggestion.function}
                    </p>
                    <p className="text-sm text-green-700">
                      <span className="font-semibold">CRISPR Strategy:</span> {suggestion.crispr}
                    </p>
                  </div>
                ))}
              </div>
            </div>

            {/* Literature Abstracts */}
            <div className="bg-white rounded-lg shadow-lg p-6">
              <div className="flex items-center gap-2 mb-4">
                <BookOpen className="text-green-600" size={24} />
                <h2 className="text-xl font-bold text-gray-800">Literature References</h2>
              </div>
              <div className="bg-gray-50 rounded p-4 text-sm text-gray-700 whitespace-pre-wrap font-mono max-h-96 overflow-y-auto">
                {results.abstracts}
              </div>
            </div>
          </div>
        )}

        {/* Info Box */}
        {!results && !loading && (
          <div className="bg-blue-50 border border-blue-200 rounded-lg p-6">
            <h3 className="font-semibold text-blue-900 mb-2">How It Works</h3>
            <ol className="list-decimal list-inside space-y-2 text-blue-800 text-sm">
              <li>Enter traits like "short trees with sweet fruit"</li>
              <li>GeneBot extracts key genetic traits from your description</li>
              <li>Searches literature database for relevant research</li>
              <li>Suggests candidate genes and potential CRISPR editing strategies</li>
            </ol>
          </div>
        )}
      </div>
    </div>
  );
};

export default AppleGeneBot;
