const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const archiver = require('archiver');

const dataDir = path.join(__dirname, '../xenium_data');

// GET /api/xenium/directory
router.get('/directory', (req, res) => {
  try {
    const structure = [];
    
    const processDirectory = (dirPath, parentPath = '') => {
      const files = [];
      const items = fs.readdirSync(dirPath);
      
      for (const item of items) {
        const itemPath = path.join(dirPath, item);
        const relativePath = path.join(parentPath, item);
        const stats = fs.statSync(itemPath);
        
        if (stats.isDirectory()) {
          files.push({
            name: item,
            path: relativePath,
            type: 'directory',
            files: processDirectory(itemPath, relativePath),
          });
        } else {
          files.push({
            name: item,
            path: relativePath,
            type: 'file',
            size: stats.size,
          });
        }
      }
      return files;
    };

    const items = fs.readdirSync(dataDir);
    for (const item of items) {
      const itemPath = path.join(dataDir, item);
      if (fs.statSync(itemPath).isDirectory()) {
        structure.push({
          name: item,
          path: item,
          type: 'directory',
          files: processDirectory(itemPath, item),
        });
      }
    }
    
    res.json(structure);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// GET /api/xenium/preview
router.get('/preview', (req, res) => {
  const { path: filePath, page = 1, per_page = 100 } = req.query;
  if (!filePath) {
    return res.status(400).json({ error: 'No file path provided' });
  }

  const fullPath = path.join(dataDir, filePath);

  if (!fs.existsSync(fullPath) || !fs.statSync(fullPath).isFile()) {
    return res.status(404).json({ error: 'File not found' });
  }

  if (filePath.endsWith('.csv') || filePath.endsWith('.csv.gz')) {
    const results = [];
    const stream = filePath.endsWith('.gz')
      ? fs.createReadStream(fullPath).pipe(require('zlib').createGunzip())
      : fs.createReadStream(fullPath);

    stream.pipe(require('csv-parser')())
      .on('data', (data) => results.push(data))
      .on('end', () => {
        const totalRows = results.length;
        const startIdx = (page - 1) * per_page;
        const endIdx = Math.min(startIdx + per_page, totalRows);
        const paginatedData = results.slice(startIdx, endIdx);

        res.json({
          columns: results.length > 0 ? Object.keys(results[0]) : [],
          data: paginatedData,
          total: totalRows,
          page: Number(page),
          per_page: Number(per_page),
        });
      });
  } else if (filePath.endsWith('.h5')) {
    const { spawn } = require('child_process');
    const pythonProcess = spawn('python', [path.join(__dirname, '../py/h5_preview.py'), fullPath]);

    let dataToSend = '';
    pythonProcess.stdout.on('data', (data) => {
      dataToSend += data.toString();
    });

    pythonProcess.on('close', (code) => {
      if (code !== 0) {
        return res.status(500).json({ error: 'Failed to preview H5 file' });
      }
      res.json(JSON.parse(dataToSend));
    });
  } else {
    res.status(400).json({ error: 'Unsupported file type' });
  }
});

// GET /api/xenium/download-folder
router.get('/download-folder', (req, res) => {
  const { path: folderPath } = req.query;
  if (!folderPath) {
    return res.status(400).json({ error: 'No folder path provided' });
  }

  const fullPath = path.join(dataDir, folderPath);

  if (!fs.existsSync(fullPath) || !fs.statSync(fullPath).isDirectory()) {
    return res.status(404).json({ error: 'Folder not found' });
  }

  const archive = archiver('zip', {
    zlib: { level: 9 }
  });

  res.attachment(`${path.basename(folderPath)}.zip`);
  archive.pipe(res);
  archive.directory(fullPath, false);
  archive.finalize();
});

// GET /api/xenium/download-all
router.get('/download-all', (req, res) => {
  const archive = archiver('zip', {
    zlib: { level: 9 }
  });

  res.attachment('all_data.zip');
  archive.pipe(res);

  const items = fs.readdirSync(dataDir);
  for (const item of items) {
    const itemPath = path.join(dataDir, item);
    if (fs.statSync(itemPath).isDirectory()) {
      archive.directory(itemPath, item);
    }
  }

  archive.finalize();
});

module.exports = router;