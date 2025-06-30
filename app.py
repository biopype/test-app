import streamlit as st
import pandas as pd
import numpy as np

# Set the page title
st.set_page_config(page_title="Test App", page_icon="🚀")

st.title("🚀 Hello from Streamlit!")
st.write("This is a simple Streamlit app to test deployment.")

# User input
name = st.text_input("Enter your name:")
if name:
    st.success(f"Welcome, {name} 👋")

# Random data chart
st.subheader("Random Data Line Chart")
data = pd.DataFrame({
    'x': np.arange(10),
    'y': np.random.randn(10)
})
st.line_chart(data.set_index('x'))