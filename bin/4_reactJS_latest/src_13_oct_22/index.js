import React from "react";
import ReactDOM from "react-dom";
import "./index.css";
import App from "./App";
//import "bootstrap/dist/css/bootstrap.css";
//import "bootstrap/dist/css/bootstrap-theme.css";
//import * as serviceWorker from "./serviceWorker";
//import App from "./spring/app1_back";
//import App from "./spring/app1";

// const container = document.getElementById('root');
// const root = ReactDOM.createRoot(container);

import {createRoot} from 'react-dom/client';
const rootElement = document.getElementById('root');
const root = createRoot(rootElement);


root.render(<App />);


//ReactDOM.render(<App />, document.getElementById("root"));

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: https://bit.ly/CRA-PWA
//serviceWorker.unregister();
