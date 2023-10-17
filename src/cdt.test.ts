import { cdt } from './cdt'
import { Serie } from './Serie'

test('Test cdt', () => {

    const vertices = [
        0.0266081599, 0.4090819884,
        0.7379528073, 0.5971111416,
        0.3705247021, 0.1465424167,
        0.6697023884, 0.6058185203,
        0.9552565538, 0.9181162736,
        0.8217644458, 0.2025187065,
        0.2787617030, 0.2862533619,
        0.9995375806, 0.1181573461,
        0.2655462782, 0.8756837537,
        0.4605185791, 0.0602445130
    ]

    const edges = [
        1, 2,
        0, 4
    ]

    // The triangles that must be created
    const expectedTriangles = Serie.create({
        array: [
            0, 2, 6,
            1, 3, 2,
            5, 2, 9,
            5, 1, 2,
            0, 9, 2,
            5, 9, 7,
            3, 6, 2,
            1, 7, 4,
            1, 5, 7,
            4, 0, 3,
            3, 1, 4,
            8, 0, 4,
            0, 6, 3
        ], itemSize: 3
    })

    const mesh = cdt({
        positions: Serie.create({ array: vertices, itemSize: 2 /* can be 3 in 3D */}),
        indices: Serie.create({ array: edges, itemSize: 2 })
    })
    
    mesh.indices.forEach((triangle: number[], i: number) => {
        const expectedTriangle = expectedTriangles.itemAt(i)
        for (let j = 0; j < 3; ++j) {
            expect(triangle[j]).toEqual(expectedTriangle[j])
        }
    })
})