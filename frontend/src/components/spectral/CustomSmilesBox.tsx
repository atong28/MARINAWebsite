import './SmilesSearchBox.css'

interface CustomSmilesBoxProps {
  value: string
  onChange: (value: string) => void
}

function CustomSmilesBox({ value, onChange }: CustomSmilesBoxProps) {
  return (
    <div className="smiles-search-box">
      <div className="smiles-search-box__header">
        <h3>Custom SMILES</h3>
        <p>Enter any SMILES string to add custom retrieval cards for the current prediction or SMILES query.</p>
      </div>
      <div className="smiles-search-box__field">
        <label htmlFor="custom-smiles-input" className="smiles-search-box__label">
          Custom SMILES String
        </label>
        <input
          id="custom-smiles-input"
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

export default CustomSmilesBox


