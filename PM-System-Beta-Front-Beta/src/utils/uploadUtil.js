import axiosInstance from './axiosInstance';
import { ElMessage } from "element-plus";  // Import ElMessage

const uploadFileInChunks = async (file, currentFileIndex, totalFiles, email, fileType, fileNumber, fun) => {
    const chunkSize = 5 * 1024 * 1024; // 每个分片的大小：5MB
    const totalChunks = Math.ceil(file.size / chunkSize);
    let scz = 0;

    const uploadChunk = async (index) => {
        if (index >= totalChunks) {
            return;
        }

        const start = index * chunkSize;
        const end = Math.min(file.size, start + chunkSize);
        const chunk = file.slice(start, end);

        const formData = new FormData();
        formData.append('file', chunk);
        formData.append('index', index);
        formData.append('totalChunks', totalChunks);
        formData.append('fileName', file.name);
        formData.append('fileSize', file.size);
        formData.append('type', fileType);
        formData.append('number', fileNumber);
        formData.append('currentFileIndex', currentFileIndex);
        formData.append('totalFiles', totalFiles);
        formData.append('email', email);
        formData.append('fun', fun); // 添加新的 fun 字段

        try {
            const response = await axiosInstance.post('/upload/upload', formData);

            // 检查响应状态和消息
            if (response.code === 200) {
                if (response.msg === "上传中") {
                    await uploadChunk(index + 1); // 递归上传下一个分片
                    scz++;
                    if (scz % 3 === 0) {
                        ElMessage.success(response.msg);
                    }
                }
                else if (response.msg === "所有文件上传成功") {
                    ElMessage.success("文件上传成功");
                }
            } else {
                ElMessage.error(response.msg); // 显示错误消息
            }
        } catch (error) {
            console.error('上传分片时出错:', error);
            ElMessage.error('上传出错，请重试！');
        }
    };

    await uploadChunk(0); // 开始上传第一个分片
};

const uploadFiles = (files, email, fileType, fileNumber, fun) => {
    if (files.length === 0) {
        ElMessage.error("没有文件上传！");
        return;
    }

    const totalFiles = files.length;

    // 按文件大小从大到小排序
    const sortedFiles = Array.from(files).sort((a, b) => b.size - a.size);

    const uploadNextFile = (index) => {
        if (index >= totalFiles) return;

        uploadFileInChunks(sortedFiles[index], index + 1, totalFiles, email, fileType, fileNumber, fun).then(() => {
            setTimeout(() => {
                uploadNextFile(index + 1); // 在上传下一个文件前等待300ms
            }, 300);
        });
    };

    uploadNextFile(0); // 开始上传第一个文件
};

export { uploadFiles };
