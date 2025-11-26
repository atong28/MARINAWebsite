import './SmilesSearchBox.css'

interface SmilesSearchBoxProps {
  value: string
  onChange: (value: string) => void
}

function SmilesSearchBox({ value, onChange }: SmilesSearchBoxProps) {
  return (
    <div className="smiles-search-box">
      <div className="smiles-search-box__header">
        <h3>SMILES Search</h3>
        <p>Enter a SMILES string to search for similar molecules.</p>
      </div>
      <div className="smiles-search-box__field">
        <label htmlFor="smiles-input" className="smiles-search-box__label">
          SMILES String
        </label>
        <input
          id="smiles-input"
          className="smiles-search-box__input"
          type="text"
          value={value}
          onChange={(e) => onChange(e.target.value)}
          placeholder="e.g. C1=CC=CC=C1"
        />
      </div>
    </div>
  )
}

export default SmilesSearchBox


