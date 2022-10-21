import React, { Component } from "react";
// import Headercon from "./components/header";
// import Footercon from "./components/footer";
import Navbar from "./components/navbar";
import logo from "./nextrap.png";
//import AppLand from "./landing_api/Landing_page";
class PlaceHolder extends Component {
  render() {
    return (
      <div id="body-wrap">
        <div id="inside">
          <div id="header">
            <img src={logo} width="980px" alt="nextrap" />
          </div>
        </div>
        <Navbar />

        <div id="footer">
          Server developed and maintained by: SBP Lab
          <br />
          &copy;Department of Biological Sciences, IISER Mohali, Sector 81, SAS
          Nagar, 140306, Punjab, India
        </div>
      </div>
    );
  }
}

export default PlaceHolder;
