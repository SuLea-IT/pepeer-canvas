// rule/fileTypeValidator.js
import i18n from "../i18n/i18n.js";
export const allowedFileExtensions = ['.tsv.gz', '.mtx.gz'];

export function validateFileType(file) {
    const fileName = file.name.toLowerCase();
    const isValidType = allowedFileExtensions.some(ext => fileName.endsWith(ext));

    if (!isValidType) {
        console.error(i18n.global.t('fileUploadConsoleError', {
            type: fileName,
            allowedTypes: allowedFileExtensions.join(', ')
        }));
    }
    return isValidType;
}
