const express = require('express');
const router = express.Router();
const path = require('path');
const fs = require('fs');

// The base directory where the binary data files are stored.
const DATA_BASE_DIR = path.join(__dirname, '..', 'data');

router.get('/', (req, res) => {
    const dataType = req.query.file; // CORRECTED: Was req.query.data
    const geneName = req.query.gene; // e.g., 'ALC'

    // --- Input Validation ---
    if (!dataType || !/^[a-zA-Z0-9_-]+$/.test(dataType)) {
        return res.status(400).json({ error: 'Invalid data type format.' });
    }
    if (!geneName || !/^[a-zA-Z0-9_.-]+$/.test(geneName)) {
        return res.status(400).json({ error: 'Invalid gene name format.' });
    }

    // --- File Path Construction ---
    // Constructs a path like: .../data/10DPA/genes/ALC.bin
    const filePath = path.join(DATA_BASE_DIR, dataType, 'genes', `${geneName}.bin`);

    console.log(`Attempting to read gene data file: ${filePath}`);

    // --- File Access and Reading ---
    fs.access(filePath, fs.constants.F_OK, (err) => {
        if (err) {
            console.error(`Gene data file not found: ${filePath}`);
            return res.status(404).json({ error: 'Gene data file not found.' });
        }

        // Use fs.readFile for a reliable, non-streamed read.
        fs.readFile(filePath, (readErr, data) => {
            if (readErr) {
                console.error(`Error reading gene data file: ${readErr}`);
                return res.status(500).json({ error: 'Failed to read gene data file.' });
            }
            
            // Send the binary data buffer directly.
            res.setHeader('Content-Type', 'application/octet-stream');
            res.send(data);
        });
    });
});

module.exports = router;
