import React, { Component } from "react";
class Headercon extends Component {
  render() {
    return (
      <div style={headerStyle}>
        <h1 style={{ fontSize: "6vw", fontFamily: "Iowan Old Style" }}>
          NEXTRaP
        </h1>
        <h4>
          <span style={boldStyle}>N</span>omenclature of{" "}
          <span style={boldStyle}>EX</span>ons in{" "}
          <span style={boldStyle}>TR</span>anscripts and{" "}
          <span style={boldStyle}>a</span>ssociated{" "}
          <span style={boldStyle}>P</span>rotein properties{" "}
        </h4>
      </div>
    );
  }
}
const headerStyle = {
  background: "Seagreen",
  fontWeight: 500,
  color: "Honeydew",
  padding: "0.5rem",
  fontFamily: "Iowan Old Style",
  marginBottom: "0px"
};
const boldStyle = {
  color: "white",
  fontWeight: "thick"
};

export default Headercon;
