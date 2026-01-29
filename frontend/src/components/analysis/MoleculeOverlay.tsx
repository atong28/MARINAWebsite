import { useEffect, useMemo, useRef } from 'react'
import './MoleculeOverlay.css'

interface MoleculeOverlayProps {
    smiles: string
    svg?: string
    selectedBits: Set<number>
}

function MoleculeOverlay({ smiles, svg, selectedBits }: MoleculeOverlayProps) {
    const containerRef = useRef<HTMLDivElement>(null)
    const svgHtml = useMemo(() => ({ __html: svg || '' }), [svg])

    useEffect(() => {
        if (!containerRef.current) {
            console.warn('[MoleculeOverlay] Container ref is null')
            return
        }

        if (!svg) {
            console.warn('[MoleculeOverlay] No SVG provided')
            return
        }

        const allOverlayGroups = containerRef.current.querySelectorAll('.bit-overlay') as NodeListOf<SVGGElement>

        if (allOverlayGroups.length === 0) {
            console.warn('[MoleculeOverlay] No overlay groups found in SVG!', {
                svgPreview: svg.substring(0, 200),
                hasBitOverlays: svg.includes('bit-overlays'),
                hasBitOverlay: svg.includes('bit-overlay-'),
            })
            return
        }

        // Hide all overlays
        allOverlayGroups.forEach((group) => {
            group.style.display = 'none'
        })

        if (selectedBits.size === 0) return

        // Show selected overlays (query by id; usually small set)
        selectedBits.forEach((bitIndex) => {
            const selector = `#bit-overlay-${bitIndex}`
            const overlayGroup = containerRef.current?.querySelector(selector) as SVGGElement | null
            if (overlayGroup) {
                overlayGroup.style.display = 'block'
            }
        })
    }, [svg, selectedBits])

    return (
        <div className="molecule-overlay-container" ref={containerRef}>
            {svg ? (
                svg.startsWith('data:image') ? (
                    <img src={svg} alt={`Molecule: ${smiles}`} className="molecule-image" />
                ) : (
                    <div dangerouslySetInnerHTML={svgHtml} className="molecule-svg" />
                )
            ) : (
                <div className="molecule-placeholder">
                    <p>No structure available</p>
                    <small>{smiles}</small>
                </div>
            )}
        </div>
    )
}

export default MoleculeOverlay
