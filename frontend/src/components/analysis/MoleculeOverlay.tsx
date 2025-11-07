import { useEffect, useRef, useCallback } from 'react'
import './MoleculeOverlay.css'

interface MoleculeOverlayProps {
  smiles: string
  svg?: string
  selectedBits: Set<number>
}

// Color palette for different bits (matches backend)
const BIT_COLORS = [
  'rgba(16,185,129,0.7)',   // Teal/green
  'rgba(59,130,246,0.7)',   // Blue
  'rgba(168,85,247,0.7)',   // Purple
  'rgba(239,68,68,0.7)',    // Red
  'rgba(245,158,11,0.7)',   // Orange
]

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
    
    console.log('[MoleculeOverlay] showBitOverlay:', {
      bitIndex,
      selector,
      found: !!overlayGroup,
      containerHasSvg: !!containerRef.current.querySelector('svg'),
    })
    
    if (overlayGroup) {
      console.log('[MoleculeOverlay] Showing overlay for bit', bitIndex, {
        currentDisplay: overlayGroup.style.display,
        id: overlayGroup.id,
        classList: Array.from(overlayGroup.classList),
      })
      
      overlayGroup.style.display = 'block'
      
      // Apply color based on bit index
      const color = BIT_COLORS[bitIndex % BIT_COLORS.length]
      if (!color) {
        console.warn('[MoleculeOverlay] No color found for bit index', bitIndex)
        return
      }
      
      const fillColor = color.replace('0.7', '0.25')
      
      // Update colors for all elements in this overlay group
      const bonds = overlayGroup.querySelectorAll('.bit-bond')
      const atoms = overlayGroup.querySelectorAll('.bit-atom')
      
      console.log('[MoleculeOverlay] Updating colors for bit', bitIndex, {
        bondsCount: bonds.length,
        atomsCount: atoms.length,
        color,
        fillColor,
      })
      
      bonds.forEach((bond) => {
        if (bond instanceof SVGLineElement) {
          bond.setAttribute('stroke', color)
        }
      })
      
      atoms.forEach((atom) => {
        if (atom instanceof SVGCircleElement) {
          atom.setAttribute('fill', fillColor)
          atom.setAttribute('stroke', color)
        }
      })
    } else {
      console.warn('[MoleculeOverlay] Overlay group not found for bit', bitIndex, {
        selector,
        allOverlayGroups: containerRef.current.querySelectorAll('.bit-overlay').length,
        allOverlayIds: Array.from(containerRef.current.querySelectorAll('.bit-overlay')).map((el) => el.id),
      })
    }
  }, [])
  
  // Update overlay visibility based on selected bits
  useEffect(() => {
    console.log('[MoleculeOverlay] useEffect triggered:', {
      hasContainer: !!containerRef.current,
      hasSvg: !!svg,
      svgLength: svg?.length || 0,
      selectedBitsCount: selectedBits.size,
      selectedBitsArray: Array.from(selectedBits),
    })
    
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
      
      console.log('[MoleculeOverlay] Overlay groups found:', {
        count: allOverlayGroups?.length || 0,
        ids: allOverlayGroups ? Array.from(allOverlayGroups).map((g) => g.id) : [],
      })
      
      if (!allOverlayGroups || allOverlayGroups.length === 0) {
        console.warn('[MoleculeOverlay] No overlay groups found in SVG!', {
          svgPreview: svg.substring(0, 200),
          hasBitOverlays: svg.includes('bit-overlays'),
          hasBitOverlay: svg.includes('bit-overlay-'),
        })
    }
        
        // Inspect SVG structure and overlay positions
        const svgElement = containerRef.current?.querySelector('svg')
        if (svgElement) {
          const viewBox = svgElement.getAttribute('viewBox')
          console.log('[MoleculeOverlay] SVG viewBox:', viewBox)
          
          // Get some atom text positions for comparison
          const textElements = svgElement.querySelectorAll('text')
          console.log(`[MoleculeOverlay] Found ${textElements.length} text elements (atoms)`)
          if (textElements.length > 0) {
            const firstText = textElements[0] as SVGTextElement
            console.log('[MoleculeOverlay] First atom text position:', {
              x: firstText.getAttribute('x'),
              y: firstText.getAttribute('y'),
              text: firstText.textContent,
            })
          }
          
          // Inspect overlay element positions
          if (allOverlayGroups && allOverlayGroups.length > 0) {
            const firstOverlay = allOverlayGroups[0]
            if (firstOverlay) {
              const lines = firstOverlay.querySelectorAll('.bit-bond')
              const circles = firstOverlay.querySelectorAll('.bit-atom')
              console.log(`[MoleculeOverlay] First overlay (${firstOverlay.id}) has ${lines.length} lines and ${circles.length} circles`)
              if (lines.length > 0) {
                const firstLine = lines[0] as SVGLineElement
                console.log('[MoleculeOverlay] First overlay line position:', {
                  x1: firstLine.getAttribute('x1'),
                  y1: firstLine.getAttribute('y1'),
                  x2: firstLine.getAttribute('x2'),
                  y2: firstLine.getAttribute('y2'),
                })
              }
              if (circles.length > 0) {
                const firstCircle = circles[0] as SVGCircleElement
                console.log('[MoleculeOverlay] First overlay circle position:', {
                  cx: firstCircle.getAttribute('cx'),
                  cy: firstCircle.getAttribute('cy'),
                  r: firstCircle.getAttribute('r'),
                })
              }
            }
          }
        }
    
      // Hide all overlays first
      allOverlayGroups?.forEach((group) => {
        group.style.display = 'none'
      })
      
      // Show overlays for selected bits
      if (selectedBits.size > 0) {
        console.log('[MoleculeOverlay] Showing overlays for selected bits:', Array.from(selectedBits))
        selectedBits.forEach((bitIndex) => {
          showBitOverlay(bitIndex)
        })
      } else {
        console.log('[MoleculeOverlay] No bits selected, all overlays hidden')
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
