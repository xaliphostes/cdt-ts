# 2D Constrained Delaunay Triangulation

Original code from [savithru-j](https://github.com/savithru-j/cdt-js) in JavaScript

This version is
1. in TypeScript
2. without any UI

# Building etc...
- `yarn` to install the dependencies
- `yarn build` to the Javascript library (in the `dist` dir)
- `yarn test` to test the lib

# Example
See also [this test](./src/cdt.test.ts) for the complete example.
```ts
import { cdt, Serie } from 'cdt-ts'

// All vertices (also containing the edge vertices)
const vertices = [
    0.026, 0.409,
    0.737, 0.597,
    0.370, 0.146,
    0.669, 0.605,
    0.955, 0.918,
    0.821, 0.202,
    0.278, 0.286,
    0.999, 0.118,
    0.265, 0.875,
    0.460, 0.060
]

// The 2 indices of each consrained edge
const edges = [
    1, 2, // first
    0, 4  // second
]

const mesh = cdt({
    positions: Serie.create({array: vertices, itemSize: 2 }),
    edges: Serie.create({ array: edges, itemSize: 2 })
})

mesh.positions.forEach( p => console.log(p))
mesh.indices.forEach( t => console.log(t))
```