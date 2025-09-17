<template>
  <div class="platform-step">
    <div class="platform-step-header">
      <div class="platform-step-title">
        <span>MapMyCells</span>
      </div>
      <div class="platform-step-steps">
        <el-steps
          class="custom-steps"
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
      </div>
    </div>
  </div>
  <div class="platform-container">
    <div class="platform-left">
      <div class="platform-step-span">
        <span class="platform-step-span-title">Step 1</span>
        <span class="platform-step-span-text"
          >Upload your gene expression data</span
        >
      </div>
      <div class="platform-left-select">
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
      <div class="platform-left-tip">
        <span>Data Usage & Privacy</span>
        <span
          >We do not use, retain, or aggregate any data uploaded to MapMyCells
          for internal purposes, nor do we publish your data. Our database
          administrators may access uploaded datasets for debugging and error
          resolution purposes. All files will be automatically deleted one week
          after upload. Please refrain from uploading any sensitive data,
          personally identifiable information, or protected health data that
          could compromise individual privacy. For more information, please
          refer to our Privacy Policy.</span
        >
      </div>
    </div>
    <div class="platform-right">
      <div class="platform-step-span">
        <span class="platform-step-span-title">Step 2</span>
        <span class="platform-step-span-text"
          >Select your data display type</span
        >
      </div>
      <div class="platform-left-select">
        <div class="platform-right-image">
          <img src="/UMAP.png" alt="" />
        </div>
        <DataSelect
          v-model="selectedValues[4]"
          :options="funOptions"
          tagText="Functionality"
          placeholderText="Select"
        />
        <!-- First DataSelect component -->
        <DataSelect
          v-model="selectedValues[2]"
          :options="options"
          tagText="Data set"
          placeholderText="Select"
        />
        <!-- New DataSelect component -->

        <DataInput
          v-model="Inputvalue"
          tagText="Enter your email address"
          placeholderText="Select"
          :validateFn="validateEmail"
        />
        <el-button @click="handleUpload">Upload</el-button>
      </div>
    </div>
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
.platform-tips {
  width: 20% !important;
}
.steps-card {
  max-width: 500px; /* 设置最大宽度 */
  margin: auto; /* 居中显示 */
}
.platform-left {
  display: flex;
  padding-left: 10%;
  flex-direction: column;
  height: 100%;
  justify-content: space-between;
}
.platform-left::after {
  content: "";
  position: absolute;
  bottom: 10%;
  left: 45%;
  /* right: 0; */
  height: 60%;
  background-color: #cccccc;
  transform: scaleX(1);
  transition: all 0.3s ease;
  width: 1px;
}
.platform-left :deep(.el-select__wrapper) {
  box-shadow: none;
  background: none;
}
.platform-right :deep(.el-select__wrapper) {
  box-shadow: none;
  background: none;
}

.platform-right :deep(.el-input__wrapper) {
  box-shadow: none;
  background: none;
}
.custom-steps :deep(.el-step__title) {
  font-size: 12px; /* 设置标题字体大小为12px */
}
.platform-step-header {
  display: flex;
}
.platform-step {
  margin-bottom: 40px;
}
.platform-step-steps {
  width: 50%;
}
.platform-step-title {
  width: 20%;
}
.platform-container {
  width: 100%;
  height: 70%;
  display: flex;
}
.platform-container > div {
  width: 40%;
}
.platform-step-title > span {
  display: flex;
  padding-left: 26%;
  font-size: 24px;
}
.platform-step-span {
  display: flex;
  font-size: 18px;
  flex-direction: column;
  margin-left: -4%;
  align-items: flex-start;
  width: 60%;
  text-align: left;
}
.platform-left-tip {
  height: 200px;
  display: flex;
  flex-direction: column;
  /* align-items: flex-start; */
  text-align: left;
  font-size: 12px;
  width: 360px;
  margin-left: -4%;
}

.platform-right {
  width: 50%;
}
.platform-right-image {
  display: flex;
  flex-direction: column-reverse;
  width: 220px;
  height: 220px;
}
.platform-right-image img {
  width: 180px;
  height: 180px;
}
.platform-tips {
  margin-top: -100px !important;
  width: 300px;
  background-color: #fff;
  border: 1px solid #ccc;
  border-radius: 5px;
  padding: 10px;
  height: calc(100vh - 96px);
  overflow-y: auto;
  display: flex;
  flex-direction: column;
  /* align-items: flex-start; */
  /* 自定义滚动条样式 */
  scrollbar-width: thin;
  /* Firefox */
  scrollbar-color: #666 transparent;
  /* Firefox */
}

.platform-tips::-webkit-scrollbar {
  /* WebKit 浏览器 */
  width: 10px;
  border-radius: 5px;
  /* 滚动条宽度 */
}

.platform-tips::-webkit-scrollbar-track {
  /* 滚动条轨道 */
  background: #f1f1f1;
}

.platform-tips::-webkit-scrollbar-thumb {
  /* 滚动条滑块 */
  background: #666;
}

.platform-tips::-webkit-scrollbar-thumb:hover {
  /* 滑块悬停状态 */
  background: #555;
}

.lj {
  padding: 0;
  margin: 0;
  line-height: 20px;
}

.conter {
  padding-left: 10%;
}

.ml {
  background: #f8f8f8;
  border-radius: 5px;
  border: 1px solid #ccc;
  padding: 10px;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
}

.title {
  color: coral;
}

.text span {
  background: #e7eaed;
  padding: 0 4px;
  border-radius: 5px;
}
</style>
