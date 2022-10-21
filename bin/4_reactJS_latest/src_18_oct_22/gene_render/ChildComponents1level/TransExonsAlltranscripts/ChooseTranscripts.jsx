import React, { Component } from "react";
import ItemsCarousel from "react-items-carousel";
// import range from "lodash/range";
// import * as Icon from 'react-bootstrap-icons'
//import { Spring } from "react-spring/renderprops";

class ChooseTranscripts extends Component {
  constructor(props) {
    super(props);
    this.state = {
      activeItemIndex: 0
    };
  }
  changeActiveItem = activeItemIndex => this.setState({ activeItemIndex });

  render() {
    const { activeItemIndex } = this.state;
    const trans = this.props.trans;
    //console.log("thispropschoose", this.props, "traNS", trans);

    let transslide = [];
    trans.forEach(tr => {
      let excount = (tr.exonsIds.match(/,/g) || []).length + 1;
      let disorder = Math.round(100 - parseFloat(tr.structured_count_disp));
      let structure = Math.round(parseFloat(tr.structured_count_ssp));
      let presence = this.props.transadded.indexOf(tr) > -1 ? true : false;
      let badgeClass = "badge ";
      let ob = (
        <div
          //className="container-fluid text-center border shadow p-2 mb-1 bg-white rounded"
          className="card border-5"
          style={{ alignItems: "center" }}
          // text-center border shadow p-2 mb-1 bg-white rounded
          key={tr.id}
        >
          <div className="flex mx-2" style={{ alignItems: "center" }}>
            <a
              href={"https://www.ncbi.nlm.nih.gov/protein/" + tr.tId}
              target="_blank"
              style={{ alignItems: "center" }}
            >
              <span
                //className="center-block rounded-pill"
                className="badge badge-pill align-center"
                style = {presence === true ? cardStyleOn : cardStyleOff}
              >
                {tr.tId}
              </span>
            </a>
          </div>
          <div
            className="text-center small"
            style={{
              padding: "0rem",
              margin: 0
            }}
          >
            <p
              style={{
                padding: "0rem",
                margin: 0
              }}
            >
              Length: {tr.length} aa {<br />}
              Exons: {excount} {<br />}Disorder: {disorder}% {<br />}
              Structure: {structure}%
            </p>
            <a href="#" style={{fontSize:"22px"}}>
              <i
                className={
                  presence === true
                    ? "bi bi-dash-circle-dotted"
                    : "bi bi-plus-circle-dotted" }
                  style = { presence === true ? buttonStyleOn : buttonStyleOff }
                
                onClick={() => this.props.handleUpdate(tr)}>
                {/* {presence === true ? "Hide" : "Show"} */}
              </i>
            </a>
          </div>
        </div>
      );
      transslide.push(ob);
    });
    return (
      <div className="container-fluid w-75 border shadow p-1 mb-2 bg-white rounded">
        <ItemsCarousel
          // Placeholder configurations
          // Carousel configurations
          numberOfCards={6}
          gutter={40}
          showSlither={false}
          slidesToScroll={3}
          firstAndLastGutter={true}
          freeScrolling={false}
          // Active item configurations
          requestToChangeActive={this.changeActiveItem}
          activeItemIndex={activeItemIndex}
          activePosition={"center"}
          chevronWidth={10}
          rightChevron={<button className="btn btn-lg btn-dark">{">"}</button>}
          leftChevron={<button className="btn btn-lg btn-dark">{"<"}</button>}
          outsideChevron={true}
        >
          {transslide}
        </ItemsCarousel>
      </div>
    );
  }
}
const cardStyleOn={
    // fontWeight: "bold",
    alignItems: "center",
    background: "#e9e9e9",
    // background: presence === true ? "peachpuff" : "mediumseagreen	",
    color: "#985e6d"
}

const cardStyleOff={
  fontWeight: "bold",
  alignItems: "center",
  background: "#e9e9e9",
  // background: presence === true ? "peachpuff" : "mediumseagreen	",
  color: "#00303f"
}

const buttonStyleOn={
  fontSize: "22px", 
  color: "#985e6d"
}

const buttonStyleOff={
  fontSize: "22px",
  color: "#00303f", 
}

export default ChooseTranscripts;

/*

*/
