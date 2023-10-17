// fmaerten (aka xaliphostes)

export class Serie {
    private itemSize_ = 2
    private array_: number[] = []

    static create({ array, itemSize }: { array: number[], itemSize: number }): Serie {
        return new Serie({ array, itemSize })
    }

    get array() {
        return this.array_
    }

    get itemSize() {
        return this.itemSize_
    }

    get count() {
        return this.array_.length / this.itemSize_
    }

    itemAt(i: number): number[] {
        const size = this.itemSize_
        const start = i * size
        const r = new Array(size).fill(0)
        for (let j = 0; j < size; ++j) {
            r[j] = this.array[start + j]
        }

        return r
    }

    forEach(callback: (item: number[], i: number, serie: Serie) => void) {
        for (let i = 0; i < this.count; ++i) {
            callback(this.itemAt(i), i, this)
        }
    }

    private constructor({ array, itemSize }: { array: number[], itemSize: number }) {
        this.array_ = array
        this.itemSize_ = itemSize
    }
}
