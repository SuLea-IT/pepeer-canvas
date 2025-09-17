# 椤圭洰鍚堝苟璁″垝锛歑enium 鏁版嵁闆嗘垚鍒?PM-System-Beta

## 1. 鍓嶇杩佺Щ

**鐩爣**: 灏?`xenium鏁版嵁` 鐨?React 缁勪欢 (`FileTree`, `DataTable`) 閲嶅啓涓?Vue 缁勪欢銆?

**姝ラ**:

1.  鍦?`PM-System-Beta-Front-Beta/src/components` 鐩綍涓嬪垱寤?`XeniumFileTree.vue` 鍜?`XeniumDataTable.vue` 鏂囦欢銆?
2.  鍙傝€?`FileTree.js` 鍜?`DataTable.js` 鐨勯€昏緫锛屼娇鐢?Vue 3 Composition API 鍜?Element Plus 缁勪欢搴撳疄鐜扮浉鍚岀殑鍔熻兘銆?
3.  鍦?`PM-System-Beta-Front-Beta/src/views` 鐩綍涓嬪垱寤轰竴涓柊鐨勮鍥剧粍浠?`XeniumData.vue`锛岀敤浜庢暣鍚?`XeniumFileTree.vue` 鍜?`XeniumDataTable.vue`銆?
4.  鍦?`PM-System-Beta-Front-Beta/src/router/index.js` 涓负 `XeniumData.vue` 娣诲姞涓€涓柊鐨勮矾鐢便€?

## 2. 鍚庣寮€鍙?

**鐩爣**: 鍦?`PM-System-Beta-Front-Beta` 鐨?Express 鍚庣瀹炵幇 `xenium鏁版嵁` 鎵€闇€鐨?API銆?

**姝ラ**:

1.  鍦?`PM-System-Beta/routes` 鐩綍涓嬪垱寤轰竴涓柊鐨勮矾鐢辨枃浠?`xenium.js`銆?
2.  鍦?`xenium.js` 涓疄鐜颁互涓嬭矾鐢憋細
    *   `GET /directory`: 鎵弿鎸囧畾鐨勬暟鎹洰褰曪紝杩斿洖鏂囦欢鍜屾枃浠跺す鐨勫眰绾х粨鏋勩€?
    *   `GET /preview`: 璇诲彇骞惰В鏋?CSV 鎴?H5 鏂囦欢锛岃繑鍥為瑙堟暟鎹€?
    *   `GET /download-folder`: 灏嗘寚瀹氭枃浠跺す鍘嬬缉涓?ZIP 鏂囦欢骞舵彁渚涗笅杞姐€?
    *   `GET /download-all`: 灏嗘墍鏈夋暟鎹枃浠跺す鍘嬬缉涓?ZIP 鏂囦欢骞舵彁渚涗笅杞姐€?
3.  鍦?`PM-System-Beta/app.js` 涓敞鍐屾柊鐨?`xenium` 璺敱銆?
4.  鍒涘缓鐩稿簲鐨?Python 鑴氭湰鏉ュ鐞?H5 鍜?CSV 鏂囦欢鐨勮В鏋愶紝骞堕€氳繃 `python-shell` 鍦?Node.js 涓皟鐢ㄣ€?

## 3. 鏁版嵁鏁村悎

**鐩爣**: 灏?`xenium鏁版嵁` 鐨勬暟鎹枃浠舵暣鍚堝埌 `PM-System-Beta-Front-Beta` 鐨勬暟鎹洰褰曚腑銆?

**姝ラ**:

1.  纭畾 `xenium鏁版嵁` 鐨勬暟鎹枃浠跺瓨鏀句綅缃€?
2.  灏嗚繖浜涙暟鎹枃浠剁Щ鍔ㄥ埌 `PM-System-Beta/data/xenium` 鐩綍銆?

## 4. UI 闆嗘垚

**鐩爣**: 灏嗘柊鐨?Vue 缁勪欢鏃犵紳闆嗘垚鍒?`PM-System-Beta-Front-Beta` 鐨勭幇鏈夊墠绔?UI 涓€?

**姝ラ**:

1.  鍦?`PM-System-Beta-Front-Beta` 鐨勫鑸爮涓坊鍔犲叆鍙ｏ紝閾炬帴鍒版柊鐨?`XeniumData` 瑙嗗浘銆?
2.  纭繚鏂伴〉闈㈢殑鏍峰紡涓庣幇鏈夌郴缁熶繚鎸佷竴鑷淬€?

## 5. 娴嬭瘯鍜岄儴缃?

**鐩爣**: 瀵归泦鎴愬悗鐨勫姛鑳借繘琛屽叏闈㈡祴璇曪紝骞跺噯澶囬儴缃层€?

**姝ラ**:

1.  娴嬭瘯鏂囦欢鍜岀洰褰曠殑鏄剧ず鏄惁姝ｇ‘銆?
2.  娴嬭瘯鏂囦欢棰勮鍔熻兘鏄惁姝ｅ父銆?
3.  娴嬭瘯鍗曚釜鏂囦欢澶瑰拰鎵€鏈夋枃浠跺す鐨勪笅杞藉姛鑳姐€?
4.  杩涜璺ㄦ祻瑙堝櫒娴嬭瘯銆?
5.  淇鍙戠幇鐨?bug銆?
6.  鍑嗗閮ㄧ讲鏂囨。銆
