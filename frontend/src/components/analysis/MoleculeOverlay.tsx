import { useEffect, useRef, useCallback } from 'react'
import './MoleculeOverlay.css'

interface MoleculeOverlayProps {
    smiles: string
    svg?: string
    selectedBits: Set<number>
}

function MoleculeOverlay({ smiles, svg, selectedBits }: MoleculeOverlayProps) {
    const containerRef = useRef<HTMLDivElement>(null)

    // Show overlay for a specific bit
    const showBitOverlay = useCallback((bitIndex: number) => {
        if (!containerRef.current) {
            console.warn('[MoleculeOverlay] showBitOverlay: containerRef.current is null')
            return
        }

        const selector = `#bit-overlay-${bitIndex}`
        const overlayGroup = containerRef.current.querySelector(selector) as SVGGElement | null

        if (overlayGroup) {
            overlayGroup.style.display = 'block'
        } else {
            console.warn('[MoleculeOverlay] Overlay group not found for bit', bitIndex, {
                selector,
                allOverlayGroups: containerRef.current.querySelectorAll('.bit-overlay').length,
                allOverlayIds: Array.from(containerRef.current.querySelectorAll('.bit-overlay')).map((el) => el.id),
            })
        }
    }, [])

    useEffect(() => {
        if (!containerRef.current) {
            console.warn('[MoleculeOverlay] Container ref is null')
            return
        }

        if (!svg) {
            console.warn('[MoleculeOverlay] No SVG provided')
            return
        }

        // Small delay to ensure SVG is rendered in DOM
        const timeoutId = setTimeout(() => {
            // Get all overlay groups
            const allOverlayGroups = containerRef.current?.querySelectorAll('.bit-overlay') as NodeListOf<SVGGElement> | undefined

            if (!allOverlayGroups || allOverlayGroups.length === 0) {
                console.warn('[MoleculeOverlay] No overlay groups found in SVG!', {
                    svgPreview: svg.substring(0, 200),
                    hasBitOverlays: svg.includes('bit-overlays'),
                    hasBitOverlay: svg.includes('bit-overlay-'),
                })
            }
            allOverlayGroups?.forEach((group) => {
                group.style.display = 'none'
            })
            if (selectedBits.size > 0) {
                selectedBits.forEach((bitIndex) => {
                    showBitOverlay(bitIndex)
                })
            }
        }, 100)

        return () => {
            clearTimeout(timeoutId)
        }
    }, [svg, selectedBits, showBitOverlay])

    return (
        <div className="molecule-overlay-container" ref={containerRef}>
            {svg ? (
                svg.startsWith('data:image') ? (
                    <img src={svg} alt={`Molecule: ${smiles}`} className="molecule-image" />
                ) : (
                    <div dangerouslySetInnerHTML={{ __html: svg }} className="molecule-svg" />
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
