const express = require('express');
const router = express.Router();
const path = require('path');
const fs = require('fs');

// The base directory where the binary data files are stored.
const DATA_BASE_DIR = path.join(__dirname, '..', 'data');

router.get('/:dataType', (req, res) => {
    const dataType = req.params.dataType;
    if (!dataType || !/^[a-zA-Z0-9_-]+$/.test(dataType)) {
        return res.status(400).json({ error: 'Invalid data type format.' });
    }

    const filePath = path.join(DATA_BASE_DIR, dataType, 'clusters.bin'); // Corrected filename to 'clusters.bin'

    fs.access(filePath, fs.constants.F_OK, (err) => {
        if (err) {
            console.error(`File not found: ${filePath}`);
            return res.status(404).json({ error: 'Data file not found.' });
        }

        fs.readFile(filePath, (readErr, data) => {
            if (readErr) {
                console.error('Error reading file:', readErr);
                return res.status(500).json({ error: 'Failed to read data file.' });
            }
            res.setHeader('Content-Type', 'application/octet-stream');
            res.send(data);
        });
    });
});

module.exports = router;
