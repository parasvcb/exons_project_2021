import React, { Fragment, useState, useEffect, useCallback } from "react";
import ProtvistaTooltip from "protvista-tooltip";
import loadWebComponent from "../utils/load-web-component";

const ProtvistaTooltipWrapper = () => {
  const [visible, setVisible] = useState(true);
  const [content, setContent] = useState("!");
  const [x, setX] = useState(200);
  const [y, setY] = useState(20);

  useEffect(() => {
    loadWebComponent("protvista-tooltip", ProtvistaTooltip);
  }, []);

  const handleToggle = useCallback(event => {
    setX(event.pageX);
    setY(event.pageY);
    setVisible(visible => !visible);
  }, []);

  const handleChangeContent = useCallback(() => {
    setContent("!".repeat(Math.ceil(Math.random() * 10)));
  });

  return (
    <Fragment>
      <button type="button" onClick={handleToggle} style={{ width: "100%" }}>
        Click to toggle tooltip targetting mouse click location
      </button>
      <button type="button" onClick={handleChangeContent}>
        Change tooltip content
      </button>
      <protvista-tooltip
        title="My tooltip"
        x={x}
        y={y}
        visible={visible ? "" : undefined}
      >
        Content of the tooltip (in <code>html</code> too{content})
      </protvista-tooltip>
      <button type="button" onClick={handleToggle} style={{ width: "100%" }}>
        Click to toggle tooltip targetting mouse click location
      </button>
    </Fragment>
  );
};

export default ProtvistaTooltipWrapper;