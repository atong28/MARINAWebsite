export const ROUTES = {
  MAIN: '/' as const,
  ANALYSIS: '/analysis' as const,
}

export type RouteKey = keyof typeof ROUTES


