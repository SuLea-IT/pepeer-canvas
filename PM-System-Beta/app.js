let createError = require('http-errors');
let express = require('express');
let path = require('path');
let cookieParser = require('cookie-parser');
let logger = require('morgan');
let cors = require('cors');
const indexRouter = require('./routes');
const clusterRouter = require('./routes/cluster');
const dataRouter = require('./routes/data');
const ipLogger = require('./util/ipLogger');
const ipRouter = require('./routes/ip');
const xeniumRouter = require('./routes/xenium');

console.log('正在初始化服务器...');

let app = express();

// 设置端口
const port = process.env.PORT || 3002;

// 简化的 CORS 配置，允许所有域名访问
const corsOptions = {
  origin: '*',  // 允许所有域名访问
  methods: ['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS'],
  allowedHeaders: ['Content-Type', 'Authorization', 'X-Requested-With'],
  credentials: false,  // 由于使用 '*'，credentials 必须设置为 false
  preflightContinue: false,
  optionsSuccessStatus: 204
};

app.use(cors(corsOptions));

// 添加预检请求处理
app.options('*', cors(corsOptions));

app.use(ipLogger);

app.use(logger('dev'));
app.use(express.json());
app.use(express.urlencoded({ extended: false }));
app.use(cookieParser());
app.use(express.static(path.join(__dirname, 'public')));
app.use('/uploads', express.static(path.join(__dirname, 'uploads')));

// API 路由
app.use('/api', function (req, res, next) {
  console.log('收到 API 请求:', req.method, req.url);
  res.header('Access-Control-Allow-Origin', '*');  // 允许所有域名访问
  res.header('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS');
  res.header('Access-Control-Allow-Headers', 'Content-Type, Authorization, Content-Length, X-Requested-With');
  next();
}, indexRouter);
app.use('/api/data', dataRouter);
app.use('/api/ip', ipRouter);
app.use('/api/xenium', xeniumRouter);

// 错误处理
app.use(function (err, req, res, next) {
  console.error('发生错误:', err.stack);
  res.status(err.status || 500);
  res.json({
    error: {
      message: err.message,
      status: err.status || 500
    }
  });
});

// 捕获未处理的异常
process.on('uncaughtException', (err) => {
  console.error('未捕获的异常:', err);
});

process.on('unhandledRejection', (reason, promise) => {
  console.error('未处理的 Promise 拒绝:', reason);
});

// 启动服务器
try {
  app.listen(port, () => {
    console.log(`服务器成功启动并运行在 http://localhost:${port}`);
    console.log('CORS 已配置为允许所有域名访问');
  });
} catch (error) {
  console.error('服务器启动失败:', error);
}

module.exports = app;
