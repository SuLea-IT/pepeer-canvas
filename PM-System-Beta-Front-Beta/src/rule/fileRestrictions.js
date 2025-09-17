export const fileTypeGenes = {
    1: {
        allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
        requiredFileNames: ["barcodes", "features", "matrix", "*", "*"],
        uploadFileCount: 5,
    },
    2: {
        allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text", ".npy"],
        requiredFileNames: [
            "barcodes",
            "features",
            "matrix",
            "barcodes_pos",
            "*",
            "*",
            "*",
        ],
        uploadFileCount: 7,
    },
    3: {
        allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
        requiredFileNames: ["barcodes", "features", "matrix", "*", "*"],
        uploadFileCount: 6,
    },
    4: {
        allowedExtensions: [".csv.gz", ".mtx.gz", ".h5", ".txt", ".text"],
        requiredFileNames: ["*", "*"],
        uploadFileCount: 4,
    },
    5: {
        allowedExtensions: [".h5ad", ".txt", ".text"],
        requiredFileNames: ["*", "*"],
        uploadFileCount: 3,
    },
};
export const fileTypeGene = {
    1: {
        allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
        requiredFileNames: ["barcodes", "features", "matrix", "*", "*"],
        uploadFileCount: 6,
    },
    2: {
        allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text", ".npy"],
        requiredFileNames: [
            "barcodes",
            "features",
            "matrix",
            "barcodes_pos",
            "*",
            "*",
            "*",
        ],
        uploadFileCount: 7,
    },
    3: {
        allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
        requiredFileNames: ["barcodes", "features", "matrix", "*", "*"],
        uploadFileCount: 6,
    },
    4: {
        allowedExtensions: [".csv.gz", ".mtx.gz", ".h5", ".txt", ".text"],
        requiredFileNames: ["*", "*"],
        uploadFileCount: 4,
    },
    5: {
        allowedExtensions: [".h5ad", ".txt", ".text"],
        requiredFileNames: ["*", "*"],
        uploadFileCount: 3,
    },
};
export const fileTypeCluster = {
    1: {
        allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
        requiredFileNames: ["barcodes", "features", "matrix", "*"],
        uploadFileCount: 4,
    },
    2: {
        allowedExtensions: [".tsv.gz", ".mtx.gz", ".npy", ".txt", ".text"],
        requiredFileNames: [
            "barcodes",
            "features",
            "matrix",
            "barcodes_pos",
            "*",
            "*",
        ],
        uploadFileCount: 6,
    },
    3: {
        allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
        requiredFileNames: ["barcodes", "features", "matrix", "*"],
        uploadFileCount: 5,
    },
    4: {
        allowedExtensions: [".csv.gz", ".mtx.gz", ".h5", ".txt", ".text"],
        requiredFileNames: ["*", "*"],
        uploadFileCount: 3,
    },
    5: {
        allowedExtensions: [".h5ad", ".txt", ".text"],
        requiredFileNames: ["*", "*"],
        uploadFileCount: 2,
    },
};