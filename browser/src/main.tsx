import { createRoot } from 'react-dom/client';
import 'molstar/build/viewer/molstar.css';
import './workspace.css';
import { ProteinInMembraneWorkspace } from './ProteinInMembraneWorkspace';

const rootElement = document.getElementById('root');
if (!rootElement) throw new Error('The workspace root is missing.');

createRoot(rootElement).render(<ProteinInMembraneWorkspace />);
