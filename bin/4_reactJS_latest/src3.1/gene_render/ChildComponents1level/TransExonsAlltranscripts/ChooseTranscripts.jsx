import React, { Component } from "react";
import ItemsCarousel from "react-items-carousel";
import range from "lodash/range";

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
          className="card"
          style={{ alignItems: "center" }}
          // text-center border shadow p-2 mb-1 bg-white rounded
          key={tr.id}
        >
          <div className="flex" style={{ alignItems: "center" }}>
            <a
              href={"https://www.ncbi.nlm.nih.gov/protein/" + tr.tId}
              target="_blank"
              style={{ alignItems: "center" }}
            >
              <span
                //className="center-block rounded-pill"
                className="badge badge-pill align-center"
                style={{
                  fontWeight: "bold",
                  alignItems: "center",
                  background:
                    presence === true ? "peachpuff" : "mediumseagreen	",
                  color: presence === true ? "black" : "white	"
                }}
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
            {/* <table
                className="table-sm font-weight-bold"
                style={{ fontSize: 10 }}
              >
                <tbody>
                  <tr>
                    <td>Length</td>
                    <td>{tr.length} aa</td>
                  </tr>
                  <tr>
                    <td>Exons</td>
                    <td>{excount}</td>
                  </tr>
                  <tr>
                    <td>Disorder</td>
                    <td>{disorder}%</td>
                  </tr>
                  <tr>
                    <td>Structure</td>
                    <td>{structure}%</td>
                  </tr>
                </tbody>
              </table> */}

            <a href="#">
              <span
                className={
                  presence === true
                    ? "badge badge-dark align-center"
                    : "badge badge-light align-center"
                }
                style={{ color: presence === true ? "white" : "black" }}
                onClick={() => this.props.handleUpdate(tr)}
              >
                {presence === true ? "Hide" : "Show"}
              </span>
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

export default ChooseTranscripts;

/*

*/
