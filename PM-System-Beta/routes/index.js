// routes/index.js

let express = require('express');
let router = express.Router();

const usersRouter = require('./users');
const uploadRouter = require('./upload');
const emailRouter = require('./email');
const jsonRouter = require('./json');
const dataRouter = require('./data');
const clusterRouter = require('./cluster'); // 引入 clusterRouter

/* GET home page. */
router.get('/', function(req, res, next) {
  res.render('index', { title: 'Express' });
});

// 配置子路由
router.use('/users', usersRouter);
router.use('/upload', uploadRouter);
router.use('/email', emailRouter);
router.use('/json', jsonRouter);
router.use('/data', dataRouter);
router.use('/cluster', clusterRouter); // 正确地将 /cluster 路径委派给 clusterRouter

module.exports = router;
