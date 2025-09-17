<template>
  <div class="content">
    <div class="banner">
      <img src="/corporate-services-overview.jpg" alt="" />
      <div class="banner-text">
        <div class="banner-title">
          <span>About LAB</span>
        </div>
      </div>
    </div>
    <div class="user-info">
      <div class="info">
        <div class="container mt-4">
          <div class="leaf">
            <div class="main-heading">
              <div class="heading">
                <h4>Message from the Director General</h4>
              </div>
            </div>
          </div>

          <div class="row mt-2">
            <div class="col-lg-4">
              <img
                src="https://www.icrisat.org/assets/template/images/director_general.jpg"
                alt=""
                class="img-fluid messagesec_image w-100"
                style="padding-top: 6px"
              />

              <h6 class="messagesec_directorname mt-2">
                <b>Dr Jacqueline d’ Arros Hughes</b>
              </h6>

              <p class="text-secondary">Director General of ICRISAT</p>
            </div>

            <div class="col-lg-8">
              <p class="message_description">
                ICRISAT stands as the premier drylands agricultural research
                institute, dedicated to uplifting smallholder farmers and
                ensuring food security in semi-arid tropics. With our
                mission-driven ethos, we are committed to elevating crop
                productivity, resilience, and sustainability vital to dryland
                communities' prosperity. This steadfast commitment profoundly
                affects millions in the world's most challenging agricultural
                regions.
              </p>

              <p>
                For over five decades, our tireless team of scientists,
                researchers, and professionals has pioneered cost-effective and
                scalable solutions empowering farming communities through
                advanced research, technology, and partnerships. From the
                world’s first pigeon pea hybrid to early maturing groundnut to
                Africa’s first biofortified pearl millet and more, our legacy
                and continued focus on key semi-arid crops, like cereals,
                pulses, and oilseeds, addresses farmers' diverse challenges –
                from adapting to climate change to bolstering soil health, water
                management, and market connections.
              </p>
            </div>
          </div>
        </div>
        <div class="container messagesec_bgcolor">
          <p>
            Leading with quiet assurance, the Institute, as a tried and tested
            leader in dry-land agri-food systems is pleased to observe its
            legacy and influence, inspiring others to emulate our impacts in the
            drylands.
          </p>
          <p>
            At ICRISAT, collaboration and knowledge exchange fuel our progress.
            Our strong and enduring partnerships with governments, research
            institutions, development agencies, and the private sector amplifies
            our impact, expediting agricultural transformation. I invite you to
            join us as a valued partner on our journey as we look to the next 50
            years of delivering tangible impacts in the drylands.
          </p>
          <img
            id="myImg"
            src="/Timeline_ICRISAT-50.jpg"
            alt=""
            class="img-fluid w-100"
          />

          <div id="myModal" class="modal" style="display: none">
            <span class="close">×</span>
            <img
              class="modal-content"
              id="img01"
              src="https://www.icrisat.org/assets/Timeline_ICRISAT-50.jpg"
            />
            <div id="caption"></div>
          </div>
        </div>
        <div class="container mt-4">
          <div class="leaf">
            <div class="main-heading">
              <div class="heading">
                <h4>Message from the Director General</h4>
              </div>
            </div>
          </div>

          <div class="row mt-2">
            <div class="col-lg-4">
              <img
                src="https://www.icrisat.org/assets/template/images/director_general.jpg"
                alt=""
                class="img-fluid messagesec_image w-100"
                style="padding-top: 6px"
              />

              <h6 class="messagesec_directorname mt-2">
                <b>Dr Jacqueline d’ Arros Hughes</b>
              </h6>

              <p class="text-secondary">Director General of ICRISAT</p>
            </div>

            <div class="col-lg-8">
              <p class="message_description">
                ICRISAT stands as the premier drylands agricultural research
                institute, dedicated to uplifting smallholder farmers and
                ensuring food security in semi-arid tropics. With our
                mission-driven ethos, we are committed to elevating crop
                productivity, resilience, and sustainability vital to dryland
                communities' prosperity. This steadfast commitment profoundly
                affects millions in the world's most challenging agricultural
                regions.
              </p>

              <p>
                For over five decades, our tireless team of scientists,
                researchers, and professionals has pioneered cost-effective and
                scalable solutions empowering farming communities through
                advanced research, technology, and partnerships. From the
                world’s first pigeon pea hybrid to early maturing groundnut to
                Africa’s first biofortified pearl millet and more, our legacy
                and continued focus on key semi-arid crops, like cereals,
                pulses, and oilseeds, addresses farmers' diverse challenges –
                from adapting to climate change to bolstering soil health, water
                management, and market connections.
              </p>
            </div>
          </div>
        </div>

        <section id="mvv">
          <div class="container mt-3">
            <div class="leaf">
              <div class="main-heading">
                <!--<hr class="side_dots me-3 ">-->

                <div class="heading">
                  <h4>Vision, Mission &amp; Values</h4>
                </div>

                <!--<hr class="side_dots_2 ms-3">-->
              </div>
            </div>
          </div>
          <div
            style="
              width: 100%;
              display: flex;
              justify-content: center;
              background: #faebd7;
            "
          >
            <div
              style="
                background-color: antiquewhite;
                padding-top: 20px;
                border: collapse;
              "
            >
              <div class="container mt-2">
                <div class="row" id="sectionBtnGroup">
                  <div class="col-lg-4 visionmission_col1">
                    <button
                      @click="currentSection = 'vision'"
                      type="button"
                      class="vision-btn visionmission_btn"
                    >
                      Vision
                    </button>
                  </div>

                  <div class="col-lg-4 visionmission_col2">
                    <button
                      @click="currentSection = 'mission'"
                      type="button"
                      class="mission-btn visionmission_btn1"
                    >
                      Mission
                    </button>
                  </div>

                  <div class="col-lg-4 visionmission_col3">
                    <button
                      @click="currentSection = 'values'"
                      type="button"
                      class="values-btn visionmission_btn1"
                    >
                      Values
                    </button>
                  </div>
                </div>

                <div class="container visionmission_bgcolor mt-5"></div>

                <div
                  class="container pb-3 mt-2"
                  id="descriptionContainer"
                  style="border: none; border-color: #f7b54f"
                >
                  <!-- Vision Section -->
                  <p
                    v-if="currentSection === 'vision'"
                    class="mt-4 visionmission_bgcolor_description visionmission_description vision"
                    style="padding-top: 22px"
                  >
                    <svg
                      class="svg-inline--fa fa-circle-check"
                      style="color: green"
                      aria-hidden="true"
                      focusable="false"
                      data-prefix="fas"
                      data-icon="circle-check"
                      role="img"
                      xmlns="http://www.w3.org/2000/svg"
                      viewBox="0 0 512 512"
                    >
                      <path
                        fill="currentColor"
                        d="M256 512c141.4 0 256-114.6 256-256S397.4 0 256 0S0 114.6 0 256S114.6 512 256 512zM369 209L241 337c-9.4 9.4-24.6 9.4-33.9 0l-64-64c-9.4-9.4-9.4-24.6 0-33.9s24.6-9.4 33.9 0l47 47L335 175c9.4-9.4 24.6-9.4 33.9 0s9.4 24.6 0 33.9z"
                      ></path>
                    </svg>
                    ICRISAT's vision is a prosperous, food-secure, and resilient
                    dryland tropics.
                  </p>

                  <!-- Mission Section -->
                  <p
                    v-if="currentSection === 'mission'"
                    class="mt-4 visionmission_bgcolor_description visionmission_description mission"
                    style="padding-top: 22px"
                  >
                    <svg
                      class="svg-inline--fa fa-circle-check"
                      style="color: green"
                      aria-hidden="true"
                      focusable="false"
                      data-prefix="fas"
                      data-icon="circle-check"
                      role="img"
                      xmlns="http://www.w3.org/2000/svg"
                      viewBox="0 0 512 512"
                    >
                      <path
                        fill="currentColor"
                        d="M256 512c141.4 0 256-114.6 256-256S397.4 0 256 0S0 114.6 0 256S114.6 512 256 512zM369 209L241 337c-9.4 9.4-24.6 9.4-33.9 0l-64-64c-9.4-9.4-9.4-24.6 0-33.9s24.6-9.4 33.9 0l47 47L335 175c9.4-9.4 24.6-9.4 33.9 0s9.4 24.6 0 33.9z"
                      ></path>
                    </svg>
                    ICRISAT’s mission is to reduce poverty, hunger, malnutrition
                    and environmental degradation in the dryland tropics.
                  </p>

                  <!-- Values Section -->
                  <p
                    v-if="currentSection === 'values'"
                    class="mt-4 visionmission_bgcolor_description visionmission_description values"
                    style="padding-top: 22px"
                  >
                    <svg
                      class="svg-inline--fa fa-circle-check"
                      style="color: green"
                      aria-hidden="true"
                      focusable="false"
                      data-prefix="fas"
                      data-icon="circle-check"
                      role="img"
                      xmlns="http://www.w3.org/2000/svg"
                      viewBox="0 0 512 512"
                    >
                      <path
                        fill="currentColor"
                        d="M256 512c141.4 0 256-114.6 256-256S397.4 0 256 0S0 114.6 0 256S114.6 512 256 512zM369 209L241 337c-9.4 9.4-24.6 9.4-33.9 0l-64-64c-9.4-9.4-9.4-24.6 0-33.9s24.6-9.4 33.9 0l47 47L335 175c9.4-9.4 24.6-9.4 33.9 0s9.4 24.6 0 33.9z"
                      ></path>
                    </svg>
                    Inclusive culture
                  </p>
                  <p
                    v-if="currentSection === 'values'"
                    class="mt-4 visionmission_bgcolor_description visionmission_description values"
                    style="padding-top: 22px"
                  >
                    <svg
                      class="svg-inline--fa fa-circle-check"
                      style="color: green"
                      aria-hidden="true"
                      focusable="false"
                      data-prefix="fas"
                      data-icon="circle-check"
                      role="img"
                      xmlns="http://www.w3.org/2000/svg"
                      viewBox="0 0 512 512"
                    >
                      <path
                        fill="currentColor"
                        d="M256 512c141.4 0 256-114.6 256-256S397.4 0 256 0S0 114.6 0 256S114.6 512 256 512zM369 209L241 337c-9.4 9.4-24.6 9.4-33.9 0l-64-64c-9.4-9.4-9.4-24.6 0-33.9s24.6-9.4 33.9 0l47 47L335 175c9.4-9.4 24.6-9.4 33.9 0s9.4 24.6 0 33.9z"
                      ></path>
                    </svg>
                    "We" not "me"
                  </p>
                </div>
              </div>
              <div class="pb-3"></div>
            </div>
          </div>
        </section>
      </div>
      <section id="of">
        <div class="container mt-2">
          <div class="leaf">
            <div class="main-heading">
              <!--<hr class="side_dots me-3 ">-->

              <div class="heading"><h4>Our Focus</h4></div>

              <!--<hr class="side_dots_2 ms-3">-->
            </div>
          </div>
        </div>

        <div class="container">
          <div class="pb-2">
            <div
              class="row"
              style="
                background-color: #d6f2ff;
                box-sizing: border-box;
                padding: 22px 10px;
              "
            >
              <div class="col-lg-11">
                <h6><b>Overcoming Poverty and Hunger</b></h6>
              </div>
              <div style="display: flex">
                <div class="col-lg-5">
                  <img
                    src="/about-poverty.jpeg"
                    alt=""
                    class="img-fluid w-100"
                    style="padding-top: 10px"
                  />
                </div>

                <div class="col-lg-7">
                  <p>
                    ICRISAT's work is dedicated to combating poverty and
                    uplifting the livelihoods of communities in the dryland
                    regions it serves. Through its research and development
                    efforts, ICRISAT focuses on improving the productivity,
                    profitability, and resilience of smallholder farmers. It
                    does so by developing and disseminating climate-smart
                    agricultural technologies, innovative farming practices, and
                    high-yielding crop varieties.
                  </p>
                  <p>
                    ICRISAT also equips farmers with the tools and knowledge to
                    enhance their agricultural productivity and generate
                    sustainable incomes. It promotes market-oriented approaches,
                    facilitating access to markets, value chains, and
                    entrepreneurship opportunities for smallholder farmers, thus
                    increasing their economic opportunities and reducing their
                    vulnerability to poverty. By empowering farmers to build
                    resilience, diversify their income sources, and improve
                    their overall economic well-being, ICRISAT's work plays a
                    vital role in overcoming poverty and creating pathways to
                    sustainable development in the dryland regions.
                  </p>
                </div>
              </div>
            </div>
          </div>

          <div class="mt-3 pb-2">
            <div
              class="row"
              style="
                background-color: #fff7e7;
                box-sizing: border-box;
                padding: 22px 10px;
              "
            >
              <div class="col-lg-11">
                <h6><b>Reducing Malnutrition</b></h6>
              </div>
              <div style="display: flex">
                <div class="col-lg-5">
                  <img
                    src="/about-poverty.jpeg"
                    alt=""
                    class="img-fluid w-100"
                    style="padding-top: 10px"
                  />
                </div>

                <div class="col-lg-7">
                  <p>
                    ICRISAT's work is instrumental in addressing the pressing
                    global challenges of hunger and malnutrition. By focusing on
                    crops that are essential to dryland communities, such as
                    chickpea, pigeonpea, groundnut, sorghum, pearl millet,
                    finger millet, small millet, and oilseeds, ICRISAT develops
                    innovative solutions to enhance crop productivity,
                    nutritional value, and accessibility. Through its research,
                    ICRISAT develops climate-resilient and nutrient-rich crop
                    varieties that can thrive in harsh conditions, providing a
                    stable food supply for vulnerable populations. Moreover,
                    ICRISAT actively promotes sustainable farming practices and
                    supports smallholder farmers in adopting improved
                    agricultural techniques to enhance food production and
                    quality.
                  </p>
                  <p>
                    By empowering communities with nutritious and diversified
                    crops, ICRISAT's efforts contribute significantly to
                    reducing hunger, improving diets, and combating
                    malnutrition.
                  </p>
                </div>
              </div>
            </div>
          </div>

          <div class="mt-3 pb-2">
            <div
              class="row"
              style="
                background-color: #f0fdff;
                box-sizing: border-box;
                padding: 22px 10px;
              "
            >
              <div class="col-lg-11">
                <h6><b>Preventing Environmental Degradation</b></h6>
              </div>
              <div style="display: flex">
                <div class="col-lg-5">
                  <img
                    src="/about-poverty.jpeg"
                    alt=""
                    class="img-fluid w-100"
                    style="padding-top: 10px"
                  />
                </div>

                <div class="col-lg-7">
                  <p>
                    ICRISAT addresses dryland environmental degradation through
                    comprehensive strategies, spanning sustainable natural
                    resource management from farm to landscape scale. Our focus
                    includes climate-smart innovations, transforming delicate
                    landscapes into resilient production systems, early warning
                    systems, and promoting best practices in landscape
                    restoration. We foster replicable solutions through
                    partnerships with both public and private sectors,
                    amplifying options for Climate-Smart Agriculture (CSA). Our
                    efforts extend to groundwater enhancement, land remediation
                    techniques, nutrient management, ecosystem services, area
                    mapping, yield estimations, and weather advisories designed
                    to fortify climate resilience.
                  </p>
                  <p>
                    An example of our innovation is the 'MRIDA' mobile gaming
                    app, designed to empower smallholder farmers in devising
                    climate-smart agricultural practices. This strategy not only
                    boosts soil carbon and yield but also cultivates resilience.
                  </p>
                </div>
              </div>
            </div>
          </div>

          <div class="mt-3 pb-2">
            <div
              class="row"
              style="
                background-color: #d6f2ff;
                box-sizing: border-box;
                padding: 22px 10px;
              "
            >
              <div class="col-lg-11">
                <h6><b>Addressing Climate Change</b></h6>
              </div>
              <div style="display: flex">
                <div class="col-lg-5">
                  <img
                    src="/about-poverty.jpeg"
                    alt=""
                    class="img-fluid w-100"
                    style="padding-top: 10px"
                  />
                </div>

                <div class="col-lg-7">
                  <p>
                    ICRISAT is at the forefront of addressing the challenges
                    posed by climate change through its innovative research and
                    development initiatives. Recognizing the significant impact
                    of climate change on agriculture, ICRISAT focuses on
                    developing climate-resilient crop varieties and implementing
                    sustainable farming practices. By leveraging advanced
                    technologies, precision agriculture, and genetic resources,
                    ICRISAT aims to enhance the adaptive capacity of smallholder
                    farmers in the face of changing climatic conditions.
                  </p>
                  <p>
                    Through its research, ICRISAT works to identify and breed
                    crops that exhibit traits such as drought tolerance, heat
                    resistance, and pest and disease resilience. Additionally,
                    ICRISAT promotes agroecological approaches, conservation
                    agriculture, and efficient water management techniques to
                    mitigate the environmental impact of farming and promote
                    sustainable resource use. ICRISAT thus contributes to
                    building agricultural resilience in the face of climate
                    change.
                  </p>
                </div>
              </div>
            </div>
          </div>
        </div>
      </section>
      <section class="pb-5" id="sp">
        <div class="container mt-2">
          <div class="leaf">
            <div class="main-heading">
              <!--<hr class="side_dots me-3 ">-->

              <div class="heading"><h4>Strategic Plan 2021-2025</h4></div>
              <br />

              <!--<hr class="side_dots_2 ms-3"><br>-->
            </div>
          </div>
        </div>

        <div class="container">
          <div class="row">
            <div class="col-lg-5">
              <img
                src="/stratergy1.png"
                alt=""
                class="w-100 strategicplan_img"
              />
            </div>

            <div class="col-lg-7">
              <p>
                ICRISAT's strategic plan for 2021-2025 embodies a bold vision to
                drive agricultural innovation and improve dryland agri-food
                systems. With a focus on enhancing productivity, resilience, and
                sustainability, the strategic plan outlines key priorities and
                pathways to address the challenges faced by dryland communities.
                It emphasizes the development of climate-smart crops, resilient
                farming practices, and efficient resource management to combat
                the impacts of climate change. The plan also underscores the
                importance of gender equality, social inclusion, and empowering
                smallholder farmers, recognizing their pivotal role in
                agricultural development.
              </p>

              <p>
                ICRISAT's strategic plan sets ambitious targets for scaling up
                research and development efforts, forging strategic
                partnerships, and leveraging digital technologies for
                transformative impact. By aligning its work with Sustainable
                Development Goals and global agendas, ICRISAT aims to catalyze
                agricultural innovation, promote sustainable food systems, and
                continue to create tangible improvements in the lives and
                livelihoods of millions of people in the drylands.
              </p>

              <p>
                <a href="#">Download a copy of the Strategic Plan</a>
              </p>
            </div>
          </div>
        </div>
      </section>
    </div>
  </div>
</template>
<script>
export default {
  data() {
    return {
      currentSection: "vision", // 默认显示 vision
    };
  },
};
</script>
<style  scoped>
* {
  padding: 0;
  margin: 0;
  font-family: "Montserrat", sans-serif;
  line-height: 1.5;
}
.banner {
  position: relative;
  width: 100%;
}
.banner img {
  height: 440px;
  width: 100%;
}
.banner-text {
  position: absolute;
  top: 0;
  right: 0;
  bottom: 0;
  left: 0;
  background-color: rgba(255, 255, 255, 0.15);
}
.banner-title {
  height: 100%;
  display: flex !important;
  justify-content: center;
  align-items: center !important;
}
.banner-title span {
  font-size: 48px;
  color: aliceblue;
  background-color: rgba(79, 79, 79, 0.3);
  padding: 3px 15px;
}
.el-main {
  padding: 0;
}
.mt-3 {
  margin: 20px 0;
}
.pb-5 {
  display: flex;
  flex-direction: column;
  align-items: center;
}
.info {
  display: flex;
  justify-content: center;
  align-items: center;
  flex-wrap: wrap;
}
.leaf {
  background-image: url(/a.png);
  background-position: center top;
  background-repeat: no-repeat;
  width: 100%;
  text-align: center;
  margin: auto;
  padding-top: 20px;
  padding-bottom: 20px;
}
.main-heading {
  display: flex;
  justify-content: center;
  align-items: center;
  width: 100%;
}
.heading {
  font-size: 30px;
  font-family: Arial, Helvetica, sans-serif;
}
.pb-2 .row {
  flex-wrap: wrap;
}
.row {
  display: -ms-flexbox;
  display: flex;
  -ms-flex-wrap: wrap;
  flex-wrap: nowrap;
}
.mt-2,
.my-2 {
  margin-top: 0.5rem !important;
}
.col-lg-4,
.col-lg-8 {
  position: relative;
  width: 100%;
  padding-right: 15px;
  padding-left: 15px;
}

.col-lg-4 {
  display: flex;
  flex-direction: column;
  align-items: flex-end;
  max-width: 33.333333%;
  text-align: left;
}
.col-lg-4 img {
  width: 50%;
}
/* .col-lg-11 {
  -ms-flex: 0 0 91.666667%;
  flex: 0 0 91.666667%;
  max-width: 91.666667%;
} */
.col-lg-5 {
  -ms-flex: 0 0 41.666667%;
  flex: 0 0 41.666667%;
  max-width: 41.666667%;
}
.col-lg-7 {
  display: flex;
  flex-direction: column;
  align-items: center;
}
.col-lg-11 {
  font-size: 24px;
}
.col-lg-7 p:last-child {
  padding: 10px 0;
}
.col-lg-7 p {
  width: 90%;
}
.col-lg-8 {
  -ms-flex: 0 0 66.666667%;
  flex: 0 0 66.666667%;
  max-width: 66.666667%;
  text-align: left;
}
.container {
  max-width: 1140px;
  text-align: left;
}

.container p {
  padding: 10px 0;
}

.w-100 {
  width: 100% !important;
}
#myImg {
  cursor: pointer;
  transition: 0.3s;
}
.visionmission_col1 {
  text-align: end;
}
[type="button"]:not(:disabled),
[type="reset"]:not(:disabled),
[type="submit"]:not(:disabled),
button:not(:disabled) {
  cursor: pointer;
}
.visionmission_btn {
  background-color: #f78a2f;
  color: #fff;
  border: none;
  width: 150px;
  padding: 15px;
  border-radius: 50px 0px 50px 0px;
  height: 60px;
}
.visionmission_btn1 {
  background-color: #fff;
  color: #f78a2f;
  border: 1px solid #f78a2f;
  width: 150px;
  padding: 15px;
  border-radius: 50px 0px 50px 0px;
  height: 60px;
}
svg:not(:host).svg-inline--fa,
svg:not(:root).svg-inline--fa {
  overflow: visible;
  box-sizing: content-box;
}
[type="button"],
[type="reset"],
[type="submit"],
button {
  -webkit-appearance: button;
}
button,
select {
  text-transform: none;
}
#mvv {
  width: 100%;
  display: flex;
  justify-content: center;
  flex-direction: column;
  align-items: center;
}
#of {
  display: flex;
  justify-content: center;
  flex-direction: column;
  align-items: center;
}
.svg-inline--fa {
  display: var(--fa-display, inline-block);
  height: 1em;
  margin-top: 5px;
  overflow: visible;
  vertical-align: -0.125em;
}
</style>