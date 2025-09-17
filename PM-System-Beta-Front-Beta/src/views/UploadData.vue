<template>
  <div class="upload-data">
    <header class="upload-data__header">
      <div class="upload-data__brand">
        <span class="upload-data__logo">MapMyCells</span>
        <p class="upload-data__tagline">
          Upload your expression matrix and watch it transform into interactive spatial insights.
        </p>
      </div>
      <el-steps
        class="upload-data__steps custom-steps"
        :active="currentStep"
        align-center
        finish-status="success"
      >
        <el-step
          v-for="(step, index) in steps"
          :key="index"
          :title="$t(`steps.${step.title}`)"
        />
      </el-steps>
    </header>

    <main class="upload-data__content">
      <section class="upload-card upload-card--primary">
        <div class="upload-card__header">
          <span class="upload-card__step">Step 1</span>
          <h2>Upload your gene expression data</h2>
          <p class="upload-card__subtitle">
            Drag &amp; drop or browse your files (tsv/gz/mtx up to 500&nbsp;MB).
          </p>
        </div>
        <div class="upload-card__uploader">
          <DataUpload
            ref="dataUploadRef"
            v-model="selectedValues[3]"
            :selected-file-type="Number(selectedValues[2])"
            :selected-fun="selectedFun"
            :file-type-cluster="fileTypeCluster"
            :file-type-gene="fileTypeGene"
            :show-file-list="showFileList"
            :upload-flag="uploadFlag"
            @click="panduan"
          />
        </div>
        <div class="upload-card__note">
          <h3>Data Usage &amp; Privacy</h3>
          <p>
            We do not use, retain, or aggregate any data uploaded to MapMyCells. Uploaded datasets are
            accessed only to resolve issues and are automatically deleted one week after upload. Avoid
            uploading sensitive or identifiable information. Refer to our privacy policy for details.
          </p>
        </div>
      </section>

      <section class="upload-card upload-card--form">
        <div class="upload-card__header">
          <span class="upload-card__step">Step 2</span>
          <h2>Select your data display type</h2>
        </div>
        <div class="upload-card__body">
          <div class="upload-card__preview">
            <img src="/UMAP.png" alt="Example spatial embedding preview" />
            <span>Preview of the cluster visualization generated after processing.</span>
          </div>
          <div class="upload-card__controls">
            <DataSelect
              class="upload-card__field"
              v-model="selectedValues[4]"
              :options="funOptions"
              tagText="Functionality"
              placeholderText="Select"
            />
            <DataSelect
              class="upload-card__field"
              v-model="selectedValues[2]"
              :options="options"
              tagText="Data set"
              placeholderText="Select"
            />
            <DataInput
              class="upload-card__field"
              v-model="Inputvalue"
              tagText="Enter your email address"
              placeholderText="name@example.com"
              :validateFn="validateEmail"
            />
            <el-button
              class="upload-card__submit"
              type="primary"
              size="large"
              @click="handleUpload"
            >
              Upload
            </el-button>
          </div>
        </div>
      </section>
    </main>

    <div class="platform-tips" style="display: none">
          <h1>使用说明</h1>
          <h2>1.配置</h2>
          <h3>a.聚类</h3>
          <p>根据自己要查询的功能准备好相应数据</p>
          <div class="ml">
            <p class="lj title">├─BTSpatial数据</p>
            <p class="lj conter">│ barcodes.tsv.gz</p>
            <p class="lj conter">│ barcodes_pos.tsv.gz</p>
            <p class="lj conter">│ features.tsv.gz</p>
            <p class="lj conter">│ matrix.mtx.gz</p>
            <p class="lj title">├─h5ad数据</p>
            <p class="lj conter">│ sc_all_adata.h5ad</p>
            <p class="lj title">├─singleCell数据</p>
            <p class="lj conter">│ barcodes.tsv.gz</p>
            <p class="lj conter">│ features.tsv.gz</p>
            <p class="lj conter">│ matrix.mtx.gz</p>
            <p class="lj title">├─singleCellSpatial数据</p>
            <p class="lj conter">│ barcodes.tsv.gz</p>
            <p class="lj conter">│ barcodes_pos.tsv.gz</p>
            <p class="lj conter">│ features.tsv.gz</p>
            <p class="lj conter">│ matrix.mtx.gz</p>
            <p class="lj title">└─Xenium数据</p>
            <p class="lj conter">│ cells.csv.gz</p>
            <p class="lj conter">│ cell_feature_matrix.h5</p>
          </div>
          <h3>b.单基因</h3>
          <div class="ml">
            <p class="lj title">├─BTSpatial数据</p>
            <p class="lj conter">│ barcodes.tsv.gz</p>
            <p class="lj conter">│ barcodes_pos.tsv.gz</p>
            <p class="lj conter">│ features.tsv.gz</p>
            <p class="lj conter">│ matrix.mtx.gz</p>
            <p class="lj title">├─h5ad数据</p>
            <p class="lj conter">│ sc_all_adata.h5ad</p>
            <p class="lj title">├─singleCell数据</p>
            <p class="lj conter">│ barcodes.tsv.gz</p>
            <p class="lj conter">│ features.tsv.gz</p>
            <p class="lj conter">│ matrix.mtx.gz</p>
            <p class="lj title">├─singleCellSpatial数据</p>
            <p class="lj conter">│ barcodes.tsv.gz</p>
            <p class="lj conter">│ S_37.npy</p>
            <p class="lj conter">│ barcodes_pos.tsv.gz</p>
            <p class="lj conter">│ features.tsv.gz</p>
            <p class="lj conter">│ matrix.mtx.gz</p>
            <p class="lj title">└─Xenium数据</p>
            <p class="lj conter">│ cells.csv.gz</p>
            <p class="lj conter">│ cell_feature_matrix.h5</p>
          </div>
          <h3>c.多基因</h3>
          <div class="ml">
            <p class="lj title">├─BTSpatial数据</p>
            <p class="lj conter">│ barcodes.tsv.gz</p>
            <p class="lj conter">│ barcodes_pos.tsv.gz</p>
            <p class="lj conter">│ features.tsv.gz</p>
            <p class="lj conter">│ matrix.mtx.gz</p>
            <p class="lj title">├─h5ad数据</p>
            <p class="lj conter">│ sc_all_adata.h5ad</p>
            <p class="lj title">├─singleCell数据</p>
            <p class="lj conter">│ barcodes.tsv.gz</p>
            <p class="lj conter">│ features.tsv.gz</p>
            <p class="lj conter">│ matrix.mtx.gz</p>
            <p class="lj title">├─singleCellSpatial数据</p>
            <p class="lj conter">│ barcodes.tsv.gz</p>
            <p class="lj conter">│ S_37.npy</p>
            <p class="lj conter">│ barcodes_pos.tsv.gz</p>
            <p class="lj conter">│ features.tsv.gz</p>
            <p class="lj conter">│ matrix.mtx.gz</p>
            <p class="lj title">└─Xenium数据</p>
            <p class="lj conter">│ cells.csv.gz</p>
            <p class="lj conter">│ cell_feature_matrix.h5</p>
          </div>
          <h2>2.具体功能说明</h2>
          <h3>a.聚类</h3>
          <p class="text">1.在<span>Functionality</span>选项卡选择“聚类”</p>
          <p class="text">2.在<span>Data set</span>选择您的数据类型</p>
          <p class="text">3.单击<span>click to upload</span>上传您的数据</p>
          <p class="text">
            4.输入邮箱后点击<span>geneid.txt</span>即可等待界面显示文件上传成功
          </p>
          <p class="text">5.系统运行结束后将会发送到您的邮箱</p>
          <h3>b.单基因映射</h3>
          <p class="text">1.在<span>Functionality</span>选项卡选择“单基因映射”</p>
          <p class="text">2.在<span>Data set</span>选择您的数据类型</p>
          <p class="text">3.单击<span>click to upload</span>上传您的数据</p>
          <p class="text">
            4.<span>geneid.txt</span>可输入多个基因<span>id</span>，结果将会为这些基因的每一个相关图例
          </p>
          <p class="text">
            5.输入邮箱后点击<span>upload</span>即可等待界面显示文件上传成功
          </p>
          <p class="text">6.系统运行结束后将会发送到您的邮箱</p>
          <h3>c.多基因映射</h3>
          <p class="text">1.在<span>Functionality</span>选项卡选择“多基因映射”</p>
          <p class="text">2.在<span>Data set</span>选择您的数据类型</p>
          <p class="text">3.单击<span>click to upload</span>上传您的数据</p>
          <p class="text">
            4.<span>geneid.txt</span>可输入多个基因<span>id</span>，结果为这些基因打分后的图例
          </p>
          <p class="text">
            5.输入邮箱后点击<span>upload</span>即可等待界面显示文件上传成功
          </p>
          <p class="text">6.系统运行结束后将会发送到您的邮箱</p>
        </div>
  </div>
</template>


<script setup>
import { ref, computed, watch } from "vue";
import { useI18n } from "vue-i18n";
import DataSelect from "../components/element/DataSelect.vue";
import DataInput from "../components/element/DataInput.vue";
import { validateEmail } from "../rule/emailValidator.js";
import DataUpload from "../components/element/DataUpload.vue";
import { uploadFiles } from "../utils/uploadUtil.js";
import { ElMessage } from "element-plus";

const currentStep = ref(1); // 当前步骤
const steps = [
  { title: "selectData", description: "" },
  { title: "uploadData", description: "" },
  { title: "configureEmail", description: "" },
];

const selectedValues = ref(["", "", "1", [], "1"]);
const Inputvalue = ref("");
const showFileList = ref(true); // To toggle file list visibility
const uploadFlag = ref(false);
let selectedFun = ref(Number(selectedValues.value[4]));
const options = [
  { value: "1", label: "Single-cell data" },
  { value: "2", label: "单细胞级别空间数据" },
  { value: "3", label: "百迈客空间转录组数据" },
  { value: "4", label: "Xenium数据" },
  { value: "5", label: "H5ad数据" },
];

const funOptions = [
  { value: "1", label: "Cluster" },
  { value: "2", label: "Single Gene Mapping" },
  { value: "3", label: "Multi-gene mapping" },
];
const fileTypeGenes = {
  1: {
    allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
    requiredFileNames: ["barcodes", "features", "matrix", "*", "*"],
    uploadFileCount: 5,
  },
  2: {
    allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text", ".npy"],
    requiredFileNames: [
      "barcodes",
      "features",
      "matrix",
      "barcodes_pos",
      "*",
      "*",
      "*",
    ],
    uploadFileCount: 7,
  },
  3: {
    allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
    requiredFileNames: ["barcodes", "features", "matrix", "*", "*"],
    uploadFileCount: 6,
  },
  4: {
    allowedExtensions: [".csv.gz", ".mtx.gz", ".h5", ".txt", ".text"],
    requiredFileNames: ["*", "*"],
    uploadFileCount: 4,
  },
  5: {
    allowedExtensions: [".h5ad", ".txt", ".text"],
    requiredFileNames: ["*", "*"],
    uploadFileCount: 3,
  },
};
const fileTypeGene = {
  1: {
    allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
    requiredFileNames: ["barcodes", "features", "matrix", "*", "*"],
    uploadFileCount: 5,
  },
  2: {
    allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text", ".npy"],
    requiredFileNames: [
      "barcodes",
      "features",
      "matrix",
      "barcodes_pos",
      "*",
      "*",
      "*",
    ],
    uploadFileCount: 6,
  },
  3: {
    allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
    requiredFileNames: ["barcodes", "features", "matrix", "*", "*"],
    uploadFileCount: 6,
  },
  4: {
    allowedExtensions: [".csv.gz", ".mtx.gz", ".h5", ".txt", ".text"],
    requiredFileNames: ["*", "*"],
    uploadFileCount: 4,
  },
  5: {
    allowedExtensions: [".h5ad", ".txt", ".text"],
    requiredFileNames: ["*", "*"],
    uploadFileCount: 3,
  },
};
const fileTypeCluster = {
  1: {
    allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
    requiredFileNames: ["barcodes", "features", "matrix", "*"],
    uploadFileCount: 4,
  },
  2: {
    allowedExtensions: [".tsv.gz", ".mtx.gz", ".npy", ".txt", ".text"],
    requiredFileNames: [
      "barcodes",
      "features",
      "matrix",
      "barcodes_pos",
      "*",
      "*",
    ],
    uploadFileCount: 6,
  },
  3: {
    allowedExtensions: [".tsv.gz", ".mtx.gz", ".txt", ".text"],
    requiredFileNames: ["barcodes", "features", "matrix", "*"],
    uploadFileCount: 5,
  },
  4: {
    allowedExtensions: [".csv.gz", ".mtx.gz", ".h5", ".txt", ".text"],
    requiredFileNames: ["*", "*"],
    uploadFileCount: 3,
  },
  5: {
    allowedExtensions: [".h5ad", ".txt", ".text"],
    requiredFileNames: ["*", "*"],
    uploadFileCount: 2,
  },
};

const fileType = ref(1);
const fileNumber = ref(4);

watch(
  () => selectedValues.value[2],
  (newValue) => {
    if (newValue) {
      let restriction;
      if (selectedFun.value == 1) {
        restriction = fileTypeCluster[Number(newValue)];
      } else if (selectedFun.value == 2) {
        restriction = fileTypeGene[Number(newValue)];
      } else if (selectedFun.value == 3) {
        restriction = fileTypeGenes[Number(selectedValues.value[2])];
      }
      if (restriction) {
        fileType.value = Number(newValue);
        fileNumber.value = restriction.uploadFileCount;

        currentStep.value = 1; // 选择了 DataSet 之后，进入步骤2：上传文件
      }
    } else {
      fileType.value = 1;
      fileNumber.value = fileTypeCluster[1].uploadFileCount;
      currentStep.value = 0; // 回到第一步
    }
  }
);
watch(
  () => selectedValues.value[4],
  (newValue) => {
    selectedFun.value = Number(newValue);

    let restriction;
    if (selectedFun.value == 1) {
      restriction = fileTypeCluster[Number(selectedValues.value[2])];
    } else if (selectedFun.value == 2) {
      restriction = fileTypeGene[Number(selectedValues.value[2])];
    } else if (selectedFun.value == 3) {
      restriction = fileTypeGenes[Number(selectedValues.value[2])];
    }

    if (restriction) {
      fileNumber.value = restriction.uploadFileCount;
    } else {
      fileNumber.value = 1;
    }
  }
);

const isFileUploaded = ref(false);

watch(Inputvalue, (newValue) => {
  if (validateEmail(newValue) && isFileUploaded.value) {
    currentStep.value = 2; // Correct email and file uploaded, move to step 3
  }
});
const panduan = () => {
  console.log(selectedValues.value[3]);
  console.log(fileNumber.value);

  if (selectedValues.value[3].length === 0) {
    console.log("判断中");

    uploadFlag.value = false;
  }
};
const handleUpload = () => {
  const email = Inputvalue.value;
  const files = selectedValues.value[3];
  // console.log(selectedValues.value[3]);

  if (files && files.length > 0) {
    if (files.length != fileNumber.value) {
      console.log(fileNumber.value);
      ElMessage.error(`请上传至少 ${fileNumber.value} 个文件。`);
      return;
    }

    const areAllFilesValid = Array.from(files).every((file) =>
      validateFile(selectedFun.value, file, fileType.value)
    );

    if (!areAllFilesValid) {
      ElMessage.error("有文件类型或名称不符合要求！");
      return;
    }

    isFileUploaded.value = true; // Mark files as uploaded

    // Correctly use validateEmail function
    if (validateEmail(email)) {
      uploadFiles(
        files,
        email,
        fileType.value,
        fileNumber.value,
        selectedFun.value
      );
      const uploadComponent = document.querySelector(
        ".upload-demo .el-upload-list"
      );
      if (uploadComponent) {
        uploadComponent.innerHTML = ""; // 清空上传的文件
        uploadFlag.value = true;
        selectedValues.value[3] = [];
      }
      console.log(selectedValues.value[3]);

      // showFileList.value = false;
    } else {
      ElMessage.error("邮箱没有填写！");
    }
  } else {
    ElMessage.error("没有文件上传！");
  }
};
const validateFile = (fun, file, fileType) => {
  let restrictions;

  if (fun == 1) {
    restrictions = fileTypeCluster[fileType];
  } else if (fun == 2) {
    restrictions = fileTypeGene[fileType];
  } else if (fun == 3) {
    restrictions = fileTypeGenes[fileType];
  }
  const fileName = file.name.toLowerCase();
  const fileExtension = `.${fileName.substring(fileName.indexOf(".") + 1)}`;

  const isValidType = restrictions.allowedExtensions.includes(fileExtension);
  const matchesRequiredName = restrictions.requiredFileNames.some(
    (name) => name === "*" || fileName.includes(name)
  );

  return isValidType && matchesRequiredName;
};
</script>



<style scoped>
.upload-data {
  --card-bg: rgba(255, 255, 255, 0.05);
  --card-border: rgba(255, 255, 255, 0.12);
  --card-shadow: 0 28px 60px rgba(15, 20, 40, 0.35);
  --text-muted: rgba(255, 255, 255, 0.7);
  min-height: calc(100vh - 120px);
  padding: clamp(24px, 4vw, 56px);
  background:
    radial-gradient(circle at 10% 20%, rgba(91, 118, 255, 0.35) 0%, rgba(17, 24, 39, 0.95) 55%),
    linear-gradient(145deg, #0f172a 0%, #111827 55%, #050911 100%);
  color: #f3f6ff;
  display: flex;
  flex-direction: column;
  gap: clamp(24px, 3vw, 40px);
}

.upload-data__header {
  display: flex;
  align-items: flex-start;
  justify-content: space-between;
  gap: clamp(16px, 3vw, 40px);
  padding: clamp(16px, 3vw, 32px);
  background: rgba(15, 23, 42, 0.6);
  border: 1px solid rgba(255, 255, 255, 0.05);
  border-radius: 20px;
  backdrop-filter: blur(18px);
  box-shadow: var(--card-shadow);
}

.upload-data__brand {
  max-width: 42ch;
  display: flex;
  flex-direction: column;
  gap: 10px;
}

.upload-data__logo {
  font-size: clamp(26px, 3vw, 36px);
  font-weight: 700;
  letter-spacing: 0.04em;
  text-transform: uppercase;
}

.upload-data__tagline {
  margin: 0;
  font-size: clamp(14px, 1.5vw, 18px);
  color: var(--text-muted);
  line-height: 1.6;
}

.upload-data__steps {
  flex: 1;
  min-width: 240px;
}

.upload-data__content {
  display: flex;
  gap: clamp(24px, 3vw, 40px);
  align-items: stretch;
  flex-wrap: wrap;
}

.upload-card {
  position: relative;
  flex: 1 1 320px;
  background: var(--card-bg);
  border: 1px solid var(--card-border);
  border-radius: 24px;
  padding: clamp(24px, 3vw, 36px);
  box-shadow: var(--card-shadow);
  display: flex;
  flex-direction: column;
  gap: 24px;
  overflow: hidden;
}

.upload-card::before {
  content: "";
  position: absolute;
  inset: 0;
  pointer-events: none;
  border-radius: inherit;
  background: linear-gradient(130deg, rgba(85, 145, 255, 0.18), transparent 70%);
  opacity: 0.7;
}

.upload-card > * {
  position: relative;
  z-index: 1;
}

.upload-card--primary::before {
  background: linear-gradient(140deg, rgba(45, 212, 191, 0.22), transparent 70%);
}

.upload-card__header {
  display: flex;
  flex-direction: column;
  gap: 10px;
}

.upload-card__step {
  font-size: 14px;
  font-weight: 600;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: rgba(148, 163, 184, 0.9);
  display: inline-flex;
  align-items: center;
  gap: 8px;
}

.upload-card__step::before {
  content: "";
  width: 28px;
  height: 2px;
  border-radius: 999px;
  background: rgba(148, 163, 184, 0.45);
}

.upload-card__header h2 {
  margin: 0;
  font-size: clamp(20px, 2.5vw, 28px);
  font-weight: 600;
  color: #ffffff;
}

.upload-card__subtitle {
  margin: 0;
  font-size: 14px;
  color: var(--text-muted);
  max-width: 38ch;
  line-height: 1.6;
}

.upload-card__uploader {
  border: 1px dashed rgba(148, 163, 184, 0.5);
  border-radius: 18px;
  padding: clamp(16px, 3vw, 24px);
  background: rgba(15, 23, 42, 0.35);
  transition: border-color 0.3s ease, background 0.3s ease;
}

.upload-card__uploader:hover {
  border-color: rgba(99, 102, 241, 0.7);
  background: rgba(15, 23, 42, 0.55);
}

.upload-card__uploader :deep(.upload-demo) {
  width: 100%;
  margin: 0;
}

.upload-card__uploader :deep(.el-upload-dragger) {
  background: transparent;
  border: none;
  padding: clamp(24px, 4vw, 40px) clamp(16px, 3vw, 32px);
  color: #e2e8f0;
}

.upload-card__uploader :deep(.el-icon--upload) {
  color: #60a5fa;
  font-size: 42px;
}

.upload-card__uploader :deep(.el-upload__text) {
  font-size: 14px;
  color: rgba(226, 232, 240, 0.85);
}

.upload-card__uploader :deep(.el-upload__tip) {
  color: rgba(125, 133, 156, 0.95);
  font-size: 12px;
}

.upload-card__note {
  background: rgba(15, 23, 42, 0.55);
  border-radius: 16px;
  padding: 18px 20px;
  border: 1px solid rgba(148, 163, 184, 0.2);
  color: var(--text-muted);
  line-height: 1.6;
  font-size: 13px;
  display: flex;
  flex-direction: column;
  gap: 8px;
}

.upload-card__note h3 {
  margin: 0;
  font-size: 14px;
  font-weight: 600;
  color: #ffffff;
  letter-spacing: 0.04em;
}

.upload-card__body {
  display: flex;
  gap: clamp(16px, 2vw, 28px);
  align-items: stretch;
  flex-wrap: wrap;
}

.upload-card__preview {
  flex: 0 0 200px;
  background: rgba(15, 23, 42, 0.55);
  border: 1px solid rgba(148, 163, 184, 0.2);
  border-radius: 18px;
  padding: 18px;
  text-align: center;
  display: flex;
  flex-direction: column;
  gap: 12px;
  justify-content: center;
  color: var(--text-muted);
  font-size: 12px;
  line-height: 1.5;
}

.upload-card__preview img {
  width: 100%;
  max-width: 180px;
  height: auto;
  border-radius: 12px;
  box-shadow: 0 15px 35px rgba(15, 23, 42, 0.4);
}

.upload-card__controls {
  flex: 1 1 240px;
  display: flex;
  flex-direction: column;
  gap: 16px;
}

.upload-card__controls :deep(.platform-left-tag) {
  width: 100%;
  margin: 0;
  gap: 10px;
}

.upload-card__controls :deep(.platform-left-tag .el-tag) {
  align-self: flex-start;
  border: none;
  background: rgba(99, 102, 241, 0.2);
  color: #c7d2fe;
  border-radius: 999px;
  padding: 4px 16px;
  font-weight: 600;
}

.upload-card__controls :deep(.platform-left-tag .el-select),
.upload-card__controls :deep(.platform-left-tag .el-input) {
  width: 100% !important;
}

.upload-card__controls :deep(.el-select__wrapper),
.upload-card__controls :deep(.el-input__wrapper) {
  background: rgba(15, 23, 42, 0.6);
  border-radius: 14px;
  border: 1px solid rgba(148, 163, 184, 0.25);
  box-shadow: none;
  transition: border-color 0.2s ease, box-shadow 0.2s ease;
}

.upload-card__controls :deep(.el-select__wrapper:hover),
.upload-card__controls :deep(.el-input__wrapper.is-focus),
.upload-card__controls :deep(.el-input__wrapper:hover) {
  border-color: rgba(99, 102, 241, 0.8);
  box-shadow: 0 0 0 2px rgba(99, 102, 241, 0.25);
}

.upload-card__controls :deep(.el-select__placeholder),
.upload-card__controls :deep(.el-input__inner) {
  color: rgba(226, 232, 240, 0.82);
}

.upload-card__controls :deep(.error-message) {
  color: #fda4af;
  font-weight: 500;
}

.upload-card__submit {
  align-self: flex-start;
  padding: 10px 28px;
  border-radius: 999px;
  box-shadow: 0 12px 30px rgba(99, 102, 241, 0.35);
  font-weight: 600;
  letter-spacing: 0.02em;
}

.upload-card__submit:hover {
  transform: translateY(-1px);
  box-shadow: 0 18px 36px rgba(99, 102, 241, 0.45);
}

.custom-steps :deep(.el-step__title) {
  color: rgba(226, 232, 240, 0.7);
  font-weight: 500;
}

.custom-steps :deep(.el-step__head.is-success) {
  border-color: #22d3ee;
  color: #22d3ee;
}

.custom-steps :deep(.el-step__title.is-success) {
  color: #22d3ee;
}

.custom-steps :deep(.el-step__icon) {
  background: rgba(15, 23, 42, 0.85);
  border-color: rgba(148, 163, 184, 0.5);
  color: rgba(226, 232, 240, 0.9);
}

@media (max-width: 1024px) {
  .upload-data__header {
    flex-direction: column;
    align-items: stretch;
  }

  .upload-data__steps {
    width: 100%;
  }
}

@media (max-width: 768px) {
  .upload-data {
    padding: 24px 18px 72px;
  }

  .upload-data__content {
    flex-direction: column;
  }

  .upload-card {
    padding: 24px;
  }

  .upload-card__body {
    flex-direction: column;
  }

  .upload-card__preview {
    flex: 1;
  }
}
</style>

