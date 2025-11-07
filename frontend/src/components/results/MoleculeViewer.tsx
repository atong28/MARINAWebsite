interface MoleculeViewerProps {
  svg?: string
  smiles: string
}

function MoleculeViewer({ svg, smiles }: MoleculeViewerProps) {
  if (svg) {
    if (svg.startsWith('data:image')) {
      // Base64 image
      return <img src={svg} alt={`Molecule: ${smiles}`} className="molecule-image" />
    } else {
      // SVG string
      return <div dangerouslySetInnerHTML={{ __html: svg }} className="molecule-svg" />
    }
  }
  
  return (
    <div className="molecule-placeholder">
      <p>No structure available</p>
      <small>{smiles}</small>
    </div>
  )
}

export default MoleculeViewer

