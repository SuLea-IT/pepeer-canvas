import axios from 'axios';

const API_BASE_URL = 'http://54.245.89.178:3002/api/xenium';

export const fetchDirectory = async () => {
    try {
        const response = await axios.get(`${API_BASE_URL}/directory`);
        return response.data;
    } catch (error) {
        console.error('Error fetching directory:', error);
        throw error;
    }
};

export const fetchFilePreview = async (path, page = 1, perPage = 100) => {
    try {
        const response = await axios.get(`${API_BASE_URL}/preview`, {
            params: {
                path,
                page,
                per_page: perPage
            }
        });
        return response.data;
    } catch (error) {
        console.error('Error fetching file preview:', error);
        throw error;
    }
};

export const downloadFolder = async (path) => {
    try {
        const response = await axios.get(`${API_BASE_URL}/download-folder`, {
            params: { path },
            responseType: 'blob'
        });

        const blob = new Blob([response.data]);
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `${path.split('/').pop()}.zip`;
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
        document.body.removeChild(a);

        return response.data;
    } catch (error) {
        console.error('Error downloading folder:', error);
        throw error;
    }
};