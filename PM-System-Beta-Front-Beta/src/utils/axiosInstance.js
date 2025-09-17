import axios from 'axios';

// 创建一个 axios 实例
const axiosInstance = axios.create({
    baseURL: 'http://54.245.89.178:3000/api', // 这里填写您的 API 基础路径
    // baseURL: 'http://127.0.0.1:3000/api', // 这里填写您的 API 基础路径
    timeout: 100000, // 请求超时时间
    headers: {
        'Content-Type': 'multipart/form-data', // 确保正确设置内容类型
    }
});


// 请求拦截器
axiosInstance.interceptors.request.use(
    (config) => {
        // 在发送请求之前做些什么
        return config;
    },
    (error) => {
        // 对请求错误做些什么
        return Promise.reject(error);
    }
);

// 响应拦截器
axiosInstance.interceptors.response.use(
    (response) => {
        // 对响应数据做点什么
        return response.data;
    },
    (error) => {
        // 对响应错误做点什么
        console.error('请求错误: ', error.response ? error.response.data : error.message);
        return Promise.reject(error);
    }
);

export default axiosInstance;
