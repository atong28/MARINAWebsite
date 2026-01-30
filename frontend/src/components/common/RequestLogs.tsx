import { useEffect, useMemo, useState } from 'react'
import './RequestLogs.css'

interface LogEntry {
  ts: number
  level: string
  logger: string
  message: string
}

interface RequestLogsProps {
  requestId?: string | null
  title?: string
}

function RequestLogs({ requestId, title = 'Request Logs' }: RequestLogsProps) {
  const [logs, setLogs] = useState<LogEntry[]>([])

  const streamUrl = useMemo(() => {
    if (!requestId) return null
    return `/api/logs/stream/${encodeURIComponent(requestId)}`
  }, [requestId])

  useEffect(() => {
    if (!streamUrl) {
      setLogs([])
      return
    }
    setLogs([])
    const source = new EventSource(streamUrl)
    source.onmessage = (event) => {
      try {
        const entry = JSON.parse(event.data) as LogEntry
        setLogs((prev) => [...prev, entry])
      } catch {
        // Ignore malformed log entries
      }
    }
    source.onerror = () => {
      source.close()
    }
    return () => {
      source.close()
    }
  }, [streamUrl])

  if (!requestId) {
    return null
  }

  return (
    <section className="request-logs">
      <div className="request-logs__header">
        <span>{title}</span>
        <span className="request-logs__id">ID: {requestId}</span>
      </div>
      <div className="request-logs__body">
        {logs.length === 0 ? (
          <div className="request-logs__empty">Waiting for logs...</div>
        ) : (
          logs.map((entry, idx) => (
            <div key={`${entry.ts}-${idx}`} className={`request-logs__line level-${entry.level.toLowerCase()}`}>
              <span className="request-logs__ts">
                {new Date(entry.ts * 1000).toLocaleTimeString()}
              </span>
              <span className="request-logs__level">{entry.level}</span>
              <span className="request-logs__message">{entry.message}</span>
            </div>
          ))
        )}
      </div>
    </section>
  )
}

export default RequestLogs
