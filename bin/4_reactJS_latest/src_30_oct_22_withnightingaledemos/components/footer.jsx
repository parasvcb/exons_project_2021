import React, { Component } from "react";
class Footercon extends Component {
  render() {
    return (
      <div style={headerStyle} className="absolute-bottom">
        <span>
          Server developed and maintained by: Paras Verma <br />{" "}
          &copy;Department of Biological Sciences, IISER Mohali, Sector 81, SAS
          Nagar, 140306, Punjab, India
        </span>
      </div>
    );
  }
}
const headerStyle = {
  background: "Seagreen",
  opacity: 0.5,
  color: "Honeydew",
  fontFamily: "font-family: Arial, Helvetica, sans-serif",
  justifyContent: "center",
  display: "flex",
  alignItems: "center",
  fontSize: "small",
  bottom: "0",
  position: "relative"
};

export default Footercon;
