const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');

// 定义数据目录的基础路径
const DATA_DIR = path.join(__dirname, '..', 'data');

/**
 * @api {get} /api/data/types Get Available Data Types
 * @apiName GetDataTypes
 * @apiGroup Data
 *
 * @apiDescription 获取 `PM-System-Beta/data` 目录下所有可用的数据类型（即子目录列表）。
 *
 * @apiSuccess {Boolean} success 请求是否成功。
 * @apiSuccess {String[]} data 包含所有数据类型名称的数组。
 *
 * @apiSuccessExample {json} Success-Response:
 *     HTTP/1.1 200 OK
 *     {
 *       "success": true,
 *       "data": ["7DPA_BMKSC", "10DPA"]
 *     }
 *
 * @apiError {Boolean} success 请求是否成功 (false)。
 * @apiError {String} message 错误信息。
 *
 * @apiErrorExample {json} Error-Response:
 *     HTTP/1.1 500 Internal Server Error
 *     {
 *       "success": false,
 *       "message": "读取数据目录时出错"
 *     }
 */
router.get('/types', (req, res) => {
  fs.readdir(DATA_DIR, { withFileTypes: true }, (err, files) => {
    if (err) {
      console.error('读取数据目录时出错:', err);
      return res.status(500).json({ success: false, message: '读取数据目录时出错' });
    }

    // 过滤出所有子目录
    const directories = files
      .filter(dirent => dirent.isDirectory())
      .map(dirent => dirent.name);

    res.json({ success: true, data: directories });
  });
});

router.post('/gene-set', async (req, res) => {
  const { file } = req.query;
  const { genes } = req.body;

  if (!file || !genes || !Array.isArray(genes) || genes.length === 0) {
    return res.status(400).json({ success: false, message: 'Missing file or genes data' });
  }

  const geneDir = path.join(DATA_DIR, file, 'genes');
  const pointData = {};

  try {
    for (const gene of genes) {
      const genePath = path.join(geneDir, `${gene}.bin`);
      if (fs.existsSync(genePath)) {
        const buffer = await fs.promises.readFile(genePath);
        const dataView = new DataView(buffer.buffer);
        const numPoints = buffer.byteLength / 12;

        for (let i = 0; i < numPoints; i++) {
          const offset = i * 12;
          const x = dataView.getFloat32(offset, true);
          const y = dataView.getFloat32(offset + 4, true);
          const value = dataView.getFloat32(offset + 8, true);
          const key = `${x},${y}`;

          if (!pointData[key]) {
            pointData[key] = { x, y, totalValue: 0, count: 0 };
          }
          pointData[key].totalValue += value;
          pointData[key].count++;
        }
      }
    }

    const numResultPoints = Object.keys(pointData).length;
    if (numResultPoints === 0) {
      return res.status(404).json({ success: false, message: 'No data found for the given genes' });
    }

    const resultBuffer = new ArrayBuffer(numResultPoints * 12);
    const resultView = new DataView(resultBuffer);
    let i = 0;

    for (const key in pointData) {
      const point = pointData[key];
      const avgValue = point.totalValue / point.count;
      const offset = i * 12;
      resultView.setFloat32(offset, point.x, true);
      resultView.setFloat32(offset + 4, point.y, true);
      resultView.setFloat32(offset + 8, avgValue, true);
      i++;
    }

    res.setHeader('Content-Type', 'application/octet-stream');
    res.send(Buffer.from(resultBuffer));

  } catch (error) {
    console.error('Error processing gene set:', error);
    res.status(500).json({ success: false, message: 'Failed to process gene set' });
  }
});

/**
 * @api {get} /api/data/:dataType/genes Get Sample Genes
 * @apiName GetSampleGenes
 * @apiGroup Data
 *
 * @apiDescription 获取指定数据类型下的所有可用基因名称列表。
 *
 * @apiParam {String} dataType 数据类型的名称 (e.g., "10DPA").
 *
 * @apiSuccess {Boolean} success 请求是否成功。
 * @apiSuccess {String[]} data 包含所有基因名称的数组。
 *
 * @apiSuccessExample {json} Success-Response:
 *     HTTP/1.1 200 OK
 *     {
 *       "success": true,
 *       "data": ["GENE1", "GENE2", "GENE3"]
 *     }
 */
router.get('/:dataType/genes', (req, res) => {
  const { dataType } = req.params;
  const { search } = req.query;

  if (!dataType) {
    return res.status(400).json({ success: false, message: 'Data type is required' });
  }

  const genesDir = path.join(DATA_DIR, dataType, 'genes');

  fs.readdir(genesDir, (err, files) => {
    if (err) {
      if (err.code === 'ENOENT') {
        return res.json({ success: true, data: [] });
      }
      console.error(`读取基因目录时出错: ${genesDir}`, err);
      return res.status(500).json({ success: false, message: '读取基因目录时出错' });
    }

    let geneNames;

    if (search) {
      // 如果有搜索查询，执行头部搜索
      const searchTerm = search.toLowerCase();
      geneNames = files
        .filter(file => file.endsWith('.bin') && file.toLowerCase().startsWith(searchTerm))
        .map(file => file.slice(0, -4));
    } else {
      // 如果没有搜索，则按文件大小排序
      const filesWithStats = files
        .filter(file => file.endsWith('.bin'))
        .map(file => {
          const filePath = path.join(genesDir, file);
          try {
            const stats = fs.statSync(filePath);
            return { name: file.slice(0, -4), size: stats.size };
          } catch (e) {
            // 如果无法获取文件状态，则忽略该文件
            return null;
          }
        })
        .filter(Boolean); // 过滤掉null值

      filesWithStats.sort((a, b) => b.size - a.size);
      geneNames = filesWithStats.map(f => f.name);
    }

    // 无论如何都只返回前10条结果
    const result = geneNames.slice(0, 10);
    res.json({ success: true, data: result });
  });
});

module.exports = router;
